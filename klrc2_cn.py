#!/usr/bin/env python3
"""
KLRC2 (NKG2C) copy number from short-read WGS depth — All of Us / GCS edition.

A self-contained port of scripts/srwgs_klrc2_depth.py + klrc2_population_freq.py.
Standard library only; no pandas, no numpy, no nkc_pangenome package.

Method (unchanged from the version validated at 222/223 against HPRC assemblies):
read counts at MAPQ >= 20 over a KLRC2 window and three nearby control loci are
converted to reads per base, and copy number is read off the ratio
target/median(controls).  Only the index and the requested slices are fetched;
the CRAM is never downloaded.

Three things differ from the 1000 Genomes version, all forced by the environment
or by measurement:

  gs:// instead of https://
      htslib has a native GCS backend (hfile_gcs) driven by GCS_OAUTH_TOKEN, and
      GCS_REQUESTER_PAYS_PROJECT when the bucket bills the reader.  Access tokens
      expire, so the token is refreshed on a timer rather than read once.

  A trimmed target window
      chr12:10431117-10434225 rather than the full KLRC2 CDS span.  Per-base depth
      in five deletion homozygotes (EAS and EUR) is exactly zero across this
      window, whereas the full CDS span retains a covered block whose left edge
      falls at chr12:10434275-10434312 in every one of them.  That retained block
      is why a deletion homozygote scores ratio*2 ~ 0.20 instead of 0 on the full
      span, and why the 1000 Genomes run needed a fitted intercept.  On the
      trimmed window copy number 0 is a true zero.  The full span is measured too,
      so the two can be compared before the trimmed window is trusted.

  Calibration fitted from the data instead of hard-coded
      There is no assembly truth set inside All of Us.  The three copy-number
      clusters are far apart and heavily populated, so the affine map from ratio
      to copy number is identified by the data alone: find the two dominant modes,
      call them 1 and 2 copies, then refine.  On the 1000 Genomes cohort this
      reproduces all 3,197 published calls exactly, without using any label.
      Fitting per sequencing batch also absorbs batch effects for free.

Usage
-----
  # 1. one sample, end to end, with timing — run this first
  klrc2_cn.py selftest --cram gs://BUCKET/PATH/NA12878.cram

  # 2. measure a cohort (resumable; --shard for dsub/Cromwell fan-out)
  klrc2_cn.py measure manifest.tsv depth.tsv --jobs 16
  klrc2_cn.py measure manifest.tsv depth.$I.tsv --shard $I/$N

  # 3. fit the calibration and emit genotypes
  klrc2_cn.py call depth.tsv calls.tsv

manifest.tsv is a TSV or CSV with a `sample` column and a `cram` (or `cram_url`)
column.  Anything else is carried through untouched.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

# --------------------------------------------------------------------------
# Loci (GRCh38).  Controls are the same three used in the validated 1000
# Genomes run, all within ~280 kb of KLRC2 on chr12p13 so that they track
# regional GC and coverage.
# --------------------------------------------------------------------------

REGIONS = [
    # name              chrom    start       end    role
    ("KLRC2_clean",   "chr12", 10431117, 10434225, "target"),
    ("KLRC2_cds",     "chr12", 10431117, 10435986, "target_full"),
    ("KLRD1",         "chr12", 10308078, 10314793, "control"),
    ("KLRK1",         "chr12", 10373114, 10388810, "control"),
    ("OLR1",          "chr12", 10159880, 10172077, "control"),
]

CONTROLS = [r[0] for r in REGIONS if r[4] == "control"]
SPAN = {r[0]: r[3] - r[2] for r in REGIONS}
LOCUS = {r[0]: f"{r[1]}:{r[2]}-{r[3]}" for r in REGIONS}

FIELDS = (["sample", "status", "elapsed_s"]
          + [f"n_{r[0]}" for r in REGIONS]
          + ["ratio_clean", "ratio_cds"])


# --------------------------------------------------------------------------
# GCS access tokens
# --------------------------------------------------------------------------

class Token:
    """A GCP access token, refreshed before it can expire mid-run.

    Tokens last about an hour.  A cohort run lasts longer than that, so the
    token is re-minted on a timer.  Each samtools call is a fresh process, so
    handing the current value to subprocess at spawn time is sufficient — there
    is no long-lived connection holding a stale credential.
    """

    CMDS = [
        ["gcloud", "auth", "application-default", "print-access-token"],
        ["gcloud", "auth", "print-access-token"],
    ]

    def __init__(self, ttl: int = 1800, static: str | None = None):
        self.ttl = ttl
        self.static = static
        self._val: str | None = None
        self._at = 0.0
        self._lock = threading.Lock()

    def get(self) -> str | None:
        if self.static:
            return self.static
        with self._lock:
            now = time.monotonic()
            if self._val is not None and now - self._at < self.ttl:
                return self._val
            last = None
            for cmd in self.CMDS:
                try:
                    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                    if r.returncode == 0 and r.stdout.strip():
                        self._val, self._at = r.stdout.strip(), now
                        return self._val
                    last = r.stderr.strip()
                except (OSError, subprocess.TimeoutExpired) as exc:
                    last = str(exc)
            raise RuntimeError(
                "could not mint a GCP access token. Run `gcloud auth login` (or "
                f"pass --token) and retry. Last error: {last}")


def build_env(token: Token | None, url: str, args) -> dict:
    env = dict(os.environ)
    if url.startswith("gs://") and token is not None:
        env["GCS_OAUTH_TOKEN"] = token.get()
        if args.user_project:
            env["GCS_REQUESTER_PAYS_PROJECT"] = args.user_project
    # Pin the CRAM reference lookup. `samtools view -c` does not decode SEQ and
    # needs no reference at all, but htslib's default REF_PATH points at the EBI
    # CRAM registry, and a silent outbound fetch from a locked-down VM hangs
    # instead of failing. Fail loudly instead, unless asked otherwise.
    if not args.allow_remote_ref:
        env.setdefault("REF_PATH", args.ref_path or os.devnull)
        if args.ref_cache:
            env["REF_CACHE"] = args.ref_cache
    return env


# --------------------------------------------------------------------------
# Measurement
# --------------------------------------------------------------------------

def count_reads(url: str, region: str, args, env: dict, cwd: str) -> int | None:
    """Reads over *region* passing the MAPQ floor, or None if the query failed.

    None and 0 are different answers and are never conflated: 0 is the expected
    result for a deletion homozygote at the target locus.

    Remote reads fail transiently under concurrency — object stores reset
    connections and time out when many workers pull indexes at once. Retried
    with backoff, because a sample dropped for a transient network error is
    indistinguishable downstream from one that could not be typed.
    """
    cmd = ["samtools", "view", "-c", "-q", str(args.mapq)]
    if args.exclude_flags:
        cmd += ["--excl-flags", args.exclude_flags]
    if args.reference:
        cmd += ["-T", args.reference]
    cmd += [url, region]

    for attempt in range(args.retries + 1):
        if attempt:
            time.sleep(min(args.backoff * 2 ** (attempt - 1), 30))
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=args.timeout, env=env, cwd=cwd)
        except (subprocess.TimeoutExpired, OSError):
            continue
        if r.returncode != 0:
            continue
        out = r.stdout.strip()
        if not out:
            continue
        try:
            return int(out)
        except ValueError:
            continue
    return None


def measure_one(sample: str, url: str, args, token: Token | None) -> dict:
    """Measure one sample. Never raises; failure is reported in `status`."""
    row = {k: "" for k in FIELDS}
    row["sample"] = sample
    t0 = time.monotonic()

    try:
        env = build_env(token, url, args)
    except RuntimeError as exc:
        row["status"] = "auth_error"
        print(f"[{sample}] {exc}", file=sys.stderr)
        return row

    # One scratch directory per sample. htslib caches the remote .crai in the
    # working directory, so all five region queries share a single index fetch —
    # and the ~1.3 MB index is removed afterwards. Left in a shared CWD this
    # accumulates about 335 GB across 250,000 samples.
    counts: dict[str, int | None] = {}
    wd = tempfile.mkdtemp(prefix=f"klrc2_{sample}_", dir=args.tmp_dir or None)
    try:
        for name, *_ in REGIONS:
            counts[name] = count_reads(url, LOCUS[name], args, env, wd)
    finally:
        shutil.rmtree(wd, ignore_errors=True)

    for name, n in counts.items():
        row[f"n_{name}"] = "" if n is None else n
    row["elapsed_s"] = f"{time.monotonic() - t0:.1f}"

    ctl = [counts[c] / SPAN[c] for c in CONTROLS if counts[c] is not None]
    if not ctl:
        row["status"] = "no_controls"
        return row
    if counts["KLRC2_clean"] is None:
        row["status"] = "target_failed"
        return row

    base = sorted(ctl)[len(ctl) // 2]
    if base <= 0:
        row["status"] = "zero_control"
        return row

    row["ratio_clean"] = f"{(counts['KLRC2_clean'] / SPAN['KLRC2_clean']) / base:.6f}"
    if counts["KLRC2_cds"] is not None:
        row["ratio_cds"] = f"{(counts['KLRC2_cds'] / SPAN['KLRC2_cds']) / base:.6f}"
    row["status"] = "ok" if len(ctl) == len(CONTROLS) else "ok_partial_controls"
    return row


def read_manifest(path: str) -> list[tuple[str, str]]:
    with open(path, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        rdr = csv.DictReader(fh, delimiter="," if head.count(",") > head.count("\t") else "\t")
        cols = {c.strip().lower(): c for c in (rdr.fieldnames or [])}
        s_col = cols.get("sample") or cols.get("person_id") or cols.get("research_id")
        u_col = next((cols[k] for k in ("cram", "cram_url", "cram_uri", "path") if k in cols), None)
        if not s_col or not u_col:
            raise SystemExit(
                f"{path}: need a sample column and a cram column; found {rdr.fieldnames}")
        out = []
        for r in rdr:
            s, u = (r[s_col] or "").strip(), (r[u_col] or "").strip()
            if s and u:
                out.append((s, u))
    return out


def cmd_measure(args) -> int:
    todo = read_manifest(args.manifest)

    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        if not 0 <= i < n:
            raise SystemExit(f"--shard {args.shard}: index must be in [0,{n})")
        todo = [t for k, t in enumerate(todo) if k % n == i]
        print(f"shard {i}/{n}: {len(todo)} samples", file=sys.stderr)

    done: set[str] = set()
    if os.path.exists(args.out) and os.path.getsize(args.out) > 0:
        with open(args.out, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("status") in ("ok", "ok_partial_controls") or args.keep_failures:
                    done.add(r["sample"])
        print(f"resuming: {len(done)} samples already measured", file=sys.stderr)
    todo = [(s, u) for s, u in todo if s not in done]
    if not todo:
        print("nothing to do", file=sys.stderr)
        return 0

    token = None if args.token == "none" else Token(ttl=args.token_ttl,
                                                    static=args.token or None)
    fresh = not done
    lock = threading.Lock()
    n_ok = 0

    with open(args.out, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                           extrasaction="ignore")
        if fresh:
            w.writeheader()
            fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = [ex.submit(measure_one, s, u, args, token) for s, u in todo]
            for k, fut in enumerate(futures, 1):
                row = fut.result()
                # Written as each sample lands, so a killed job keeps its work.
                with lock:
                    w.writerow(row)
                    fh.flush()
                    if row["status"].startswith("ok"):
                        n_ok += 1
                    if k % args.report_every == 0 or k == len(todo):
                        print(f"  {k}/{len(todo)}  ok={n_ok}", flush=True, file=sys.stderr)

    print(f"wrote {args.out}: {n_ok}/{len(todo)} measured", file=sys.stderr)
    return 0 if n_ok else 1


# --------------------------------------------------------------------------
# Calibration and calling
# --------------------------------------------------------------------------

def histogram_modes(xs: list[float], bandwidth: float = 0.08,
                    min_separation: float = 0.4, nbins: int = 400):
    """Modes of *xs*, strongest first.

    The histogram is convolved with a Gaussian whose width is fixed in data
    units rather than in bins, so the result does not change with sample size:
    a cohort of 223 and a cohort of 250,000 are smoothed identically.  Binning
    first keeps the cost independent of n.

    Peaks closer together than *min_separation* are merged, since a single
    copy-number cluster can otherwise register as two adjacent bumps.
    """
    lo, hi = min(xs), max(xs)
    if hi <= lo:
        return []
    width = (hi - lo) / nbins
    h = [0.0] * nbins
    for v in xs:
        h[min(max(int((v - lo) / width), 0), nbins - 1)] += 1

    sigma = max(bandwidth / width, 1.0)
    half = int(3 * sigma)
    kernel = [2.718281828459045 ** (-0.5 * (d / sigma) ** 2)
              for d in range(-half, half + 1)]
    ktot = sum(kernel)
    sm = [0.0] * nbins
    for j in range(nbins):
        acc = 0.0
        for d, kv in enumerate(kernel, start=-half):
            i = j + d
            if 0 <= i < nbins:
                acc += h[i] * kv
        sm[j] = acc / ktot

    peaks = [(sm[j], lo + (j + 0.5) * width)
             for j in range(1, nbins - 1)
             if sm[j] > sm[j - 1] and sm[j] >= sm[j + 1] and sm[j] > 0]
    peaks.sort(reverse=True)

    kept: list[tuple[float, float]] = []
    for strength, pos in peaks:
        if all(abs(pos - p) >= min_separation for _, p in kept):
            kept.append((strength, pos))
    return kept


def calibrate(ratios: list[float], iters: int = 6) -> tuple[float, float]:
    """Fit ratio*2 = a + b*CN with no truth labels.

    The two most populated clusters are one and two copies (a deletion allele
    frequency near 0.2 puts roughly 64% of samples at two copies and 32% at
    one), which fixes the affine map.  It is then refined by alternating
    assignment and weighted least squares over the cluster means.
    """
    xs = [r * 2 for r in ratios]
    if len(xs) < 50:
        raise SystemExit(f"calibration needs >= 50 measured samples, got {len(xs)}")
    peaks = histogram_modes(xs)
    if len(peaks) < 2:
        raise SystemExit("calibration failed: fewer than two modes in the ratio "
                         "distribution — inspect the histogram (`call --hist`)")
    m1, m2 = sorted(p[1] for p in peaks[:2])
    b = m2 - m1
    if not 0.5 < b < 1.5:
        raise SystemExit(f"calibration failed: dominant modes are {b:.3f} apart, "
                         "expected ~1 copy — inspect the histogram (`call --hist`)")
    a = m1 - b

    for _ in range(iters):
        groups: dict[int, list[float]] = {}
        for x in xs:
            k = round((x - a) / b)
            if 0 <= k <= 4:
                groups.setdefault(k, []).append(x)
        floor = max(3, int(0.002 * len(xs)))
        pts = [(k, sum(g) / len(g), len(g)) for k, g in groups.items() if len(g) >= floor]
        if len(pts) < 2:
            break
        tw = sum(w for _, _, w in pts)
        mk = sum(w * k for k, _, w in pts) / tw
        mx = sum(w * m for _, m, w in pts) / tw
        den = sum(w * (k - mk) ** 2 for k, _, w in pts)
        if den == 0:
            break
        b = sum(w * (k - mk) * (m - mx) for k, m, w in pts) / den
        a = mx - b * mk
    return a, b


def cmd_call(args) -> int:
    with open(args.depth, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]

    col = "ratio_cds" if args.window == "cds" else "ratio_clean"
    # A resumed run re-measures samples that previously failed, so the file can
    # hold an old failed row and a new successful one for the same sample.
    by_sample = {}
    for r in rows:
        if r.get("status", "").startswith("ok") and r.get(col):
            by_sample[r["sample"]] = r
    usable = list(by_sample.values())
    print(f"{len(usable)}/{len({r['sample'] for r in rows})} samples with a usable {col}",
          file=sys.stderr)
    if not usable:
        return 1
    ratios = [float(r[col]) for r in usable]

    if args.calibration == "none":
        a, b = 0.0, 1.0
    elif args.calibration == "fit":
        a, b = calibrate(ratios)
    else:
        a, b = (float(v) for v in args.calibration.split(","))
    print(f"calibration: ratio*2 = {a:.4f} + {b:.4f} * CN", file=sys.stderr)
    if abs(a) > 0.35 or not 0.6 < b < 1.4:
        print("WARNING: calibration is far from the expected a~0, b~1 — check the "
              "histogram before trusting these calls", file=sys.stderr)

    if args.hist:
        for strength, pos in sorted(histogram_modes([r * 2 for r in ratios])[:6],
                                    key=lambda p: p[1]):
            print(f"  mode at ratio*2 = {pos:6.3f}   (CN {(pos - a) / b:5.2f})",
                  file=sys.stderr)

    out_fields = ["sample", "ratio", "cn_continuous", "cn", "dist", "call"]
    dist_by_cn: dict[int, list[int]] = {}
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for r, ratio in zip(usable, ratios):
            cont = (ratio * 2 - a) / b
            cn = min(max(round(cont), 0), 6)
            dist = abs(cont - cn)
            ok = dist <= args.max_dist
            dist_by_cn.setdefault(cn, []).append(int(not ok))
            w.writerow({"sample": r["sample"], "ratio": f"{ratio:.6f}",
                        "cn_continuous": f"{cont:.4f}", "cn": cn,
                        "dist": f"{dist:.4f}", "call": cn if ok else "ambiguous"})

    print(f"\n{'CN':>3} {'n':>7} {'ambiguous':>10} {'rate':>7}", file=sys.stderr)
    for cn in sorted(dist_by_cn):
        v = dist_by_cn[cn]
        print(f"{cn:>3} {len(v):>7} {sum(v):>10} {100 * sum(v) / len(v):>6.1f}%",
              file=sys.stderr)
    # A gate that fires unevenly across copy-number classes biases the allele
    # frequency. Report per class so that is visible rather than assumed away.
    rates = [sum(v) / len(v) for cn, v in dist_by_cn.items() if len(v) >= 20]
    if rates and max(rates) - min(rates) > 0.05:
        print("\nWARNING: the ambiguity gate fires at different rates across copy-number\n"
              "classes. Filtering on it will bias the deletion allele frequency; compute\n"
              "frequencies on the unfiltered calls and use the gate for QC only.",
              file=sys.stderr)

    called = sum(len(v) - sum(v) for v in dist_by_cn.values())
    n0 = len(dist_by_cn.get(0, [])); n1 = len(dist_by_cn.get(1, []))
    n = sum(len(v) for v in dist_by_cn.values())
    print(f"\nwrote {args.out}: {called}/{n} unambiguous", file=sys.stderr)
    print(f"deletion allele frequency (all calls, ungated): "
          f"{(2 * n0 + n1) / (2 * n):.4f}", file=sys.stderr)
    return 0


# --------------------------------------------------------------------------
# Selftest
# --------------------------------------------------------------------------

def cmd_selftest(args) -> int:
    if shutil.which("samtools") is None:
        print("FAIL: samtools is not on PATH", file=sys.stderr)
        return 1
    v = subprocess.run(["samtools", "--version"], capture_output=True, text=True)
    feat = next((l for l in v.stdout.splitlines() if "GCS=" in l), "")
    print(v.stdout.splitlines()[0])
    print(f"htslib features: {feat.strip() or 'unknown'}")
    if args.cram.startswith("gs://") and "GCS=yes" not in feat:
        print("FAIL: this samtools was built without GCS support, so gs:// URLs "
              "will not open. Install a samtools with GCS=yes.", file=sys.stderr)
        return 1

    token = None if args.token == "none" else Token(static=args.token or None)
    if args.cram.startswith("gs://") and token is not None:
        try:
            token.get()
            print("access token: ok")
        except RuntimeError as exc:
            print(f"FAIL: {exc}", file=sys.stderr)
            return 1

    t0 = time.monotonic()
    row = measure_one("selftest", args.cram, args, token)
    dt = time.monotonic() - t0

    print(f"\nstatus: {row['status']}   elapsed: {dt:.1f}s")
    for name, *_ in REGIONS:
        n = row[f"n_{name}"]
        d = f"{int(n) / SPAN[name]:.4f}" if n != "" else "  --  "
        print(f"  {name:14s} {LOCUS[name]:26s} reads={str(n):>8}  reads/bp={d}")
    if row["ratio_clean"]:
        rc = float(row["ratio_clean"])
        print(f"\nratio (trimmed window) = {rc:.4f}  -> {rc * 2:.2f} copies uncalibrated")
        print(f"ratio (full CDS span)  = {row['ratio_cds'] or 'n/a'}")
        print(f"\nprojected cohort throughput at --jobs {args.jobs}: "
              f"{args.jobs * 60 / max(dt, 1e-9):.0f} samples/min")
        return 0
    print("\nFAIL: no ratio produced — see status above", file=sys.stderr)
    return 1


# --------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def shared(sp):
        sp.add_argument("--mapq", type=int, default=20,
                        help="MAPQ floor; 20 excludes reads multi-mapping across "
                             "the KLRC1/KLRC2/KLRC3 paralogues [%(default)s]")
        sp.add_argument("--exclude-flags", default="",
                        help="samtools --excl-flags value, e.g. 0xF04 to drop "
                             "duplicates and supplementary reads (default: none, "
                             "matching the validated run)")
        sp.add_argument("--timeout", type=int, default=300,
                        help="per-region query timeout, seconds [%(default)s]")
        sp.add_argument("--retries", type=int, default=3,
                        help="retries per region query; remote reads fail "
                             "transiently under concurrency [%(default)s]")
        sp.add_argument("--backoff", type=float, default=2.0,
                        help="initial retry backoff, seconds [%(default)s]")
        sp.add_argument("--reference", default="",
                        help="reference FASTA; not needed for read counting, "
                             "accepted for unusual CRAMs")
        sp.add_argument("--ref-path", default="",
                        help="value for REF_PATH (default: block remote lookups)")
        sp.add_argument("--ref-cache", default="",
                        help="value for REF_CACHE")
        sp.add_argument("--allow-remote-ref", action="store_true",
                        help="permit htslib's default reference lookup, which "
                             "reaches out to the EBI CRAM registry")
        sp.add_argument("--user-project", default="",
                        help="billing project for requester-pays buckets "
                             "(sets GCS_REQUESTER_PAYS_PROJECT)")
        sp.add_argument("--token", default="",
                        help="access token to use verbatim, or 'none' to skip "
                             "token handling (default: mint via gcloud)")
        sp.add_argument("--token-ttl", type=int, default=1800,
                        help="seconds before the token is re-minted [%(default)s]")
        sp.add_argument("--tmp-dir", default="",
                        help="parent for per-sample scratch directories")
        sp.add_argument("--jobs", type=int, default=8,
                        help="concurrent samples [%(default)s]")

    m = sub.add_parser("measure", help="measure read depth for a cohort")
    m.add_argument("manifest")
    m.add_argument("out")
    m.add_argument("--shard", default="",
                   help="process shard i of n, as i/n, for batch fan-out")
    m.add_argument("--keep-failures", action="store_true",
                   help="on resume, do not retry samples that previously failed")
    m.add_argument("--report-every", type=int, default=50)
    shared(m)
    m.set_defaults(func=cmd_measure)

    c = sub.add_parser("call", help="fit the calibration and assign copy number")
    c.add_argument("depth")
    c.add_argument("out")
    c.add_argument("--window", choices=["clean", "cds"], default="clean",
                   help="which target window to call on [%(default)s]")
    c.add_argument("--calibration", default="fit",
                   help="'fit' (from the data), 'none', or literal 'a,b' "
                        "[%(default)s]")
    c.add_argument("--max-dist", type=float, default=0.18,
                   help="ambiguity gate on |CN - nearest integer|; for QC only, "
                        "not for frequency estimation [%(default)s]")
    c.add_argument("--hist", action="store_true", help="print the fitted modes")
    c.set_defaults(func=cmd_call)

    s = sub.add_parser("selftest", help="measure one sample and report timing")
    s.add_argument("--cram", required=True)
    shared(s)
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
