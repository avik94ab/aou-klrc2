#!/usr/bin/env python3
"""
Stage the NKC region out of All of Us srWGS CRAMs — GCS edition.

A port of scripts/srwgs_extract_reads.sh to gs:// paths and the Researcher
Workbench.  Standard library only; one file, so it can be copied into a
notebook or a dsub image on its own.

Capture is a ranged fetch, exactly as on the 1000 Genomes CRAMs.  The
alignments are already on GRCh38, so the two NKC intervals come out by
coordinate: ~35 MB moves per sample against a ~15 GB CRAM, and nothing close to
a whole CRAM is ever downloaded.

Per sample this writes the same three files the validated genotyper consumes,
produced by the same commands so that calls are comparable across cohorts:

    {sample}.nkc.bam      the raw slice
    {sample}.nkc.bam.bai  index; klrc2_cn.py can then measure copy number from
                          the local slice instead of re-querying the CRAM, as
                          all five of its loci fall inside the core interval
    {sample}.nkc.fq.gz    reads, names suffixed /1 /2

Three things differ from the Wynton version, all forced by the environment.

  A reference is mandatory here.
      klrc2_cn.py gets away with pinning REF_PATH at /dev/null because
      `samtools view -c` never decodes SEQ.  Writing a BAM does decode SEQ, and
      CRAM decoding needs the exact sequence the CRAM was compressed against.
      `nkc_extract.py ref` builds that cache: it reads the M5 of chr12 out of
      the CRAM header, byte-range-fetches that one contig out of a public
      GRCh38 FASTA (134 MB moved, not the 3.1 GB file), and stores it under the
      M5 htslib will look for.  Nothing outbound is permitted afterwards, so a
      reference mismatch fails loudly at cache-build time — naming the
      reference the CRAM header itself points at — instead of silently decoding
      against the wrong sequence or hanging on the EBI registry.

  gs:// instead of https://
      htslib's native GCS backend (GCS=yes in `samtools --version`) is driven
      by GCS_OAUTH_TOKEN, plus GCS_REQUESTER_PAYS_PROJECT when the bucket bills
      the reader.  Tokens expire in about an hour and a cohort run lasts
      longer, so the token is re-minted on a timer rather than read once.

  Destinations may be buckets.
      A Workbench VM or a dsub task has a small local disk and 80 MB of output
      per sample.  A gs:// destination stages each sample locally and uploads
      it, and resume is answered from one bucket listing rather than a stat per
      sample.

Usage
-----
  # 0. one-time: reference cache, keyed by the M5 in your own CRAM's header
  nkc_extract.py ref --from-cram gs://BUCKET/PATH/wgs_1000000.cram

  # 1. one sample end to end, with timing — always run this first
  nkc_extract.py selftest --cram gs://BUCKET/PATH/wgs_1000000.cram --keep reads/

  # 2. the sample you have in hand, straight from the manifest line
  nkc_extract.py extract \
      '1000000,gs://BUCKET/PATH/wgs_1000000.cram,gs://BUCKET/PATH/wgs_1000000.cram.crai' \
      reads/

  # 3. a cohort; resumable, and --shard fans out across dsub tasks
  nkc_extract.py extract manifest.tsv gs://WORKSPACE_BUCKET/nkc/reads --jobs 8
  nkc_extract.py extract manifest.tsv reads/ --shard $I/$N --jobs 8

  # 4. confirm what landed is complete and readable
  nkc_extract.py verify reads/ --manifest manifest.tsv

manifest.tsv is a TSV or CSV with a sample column, a cram column and,
optionally, a crai column; a headerless file of `sample,cram,crai` lines works
too.  Extracted reads are individual-level data: keep the destination inside
the workspace bucket.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

# --------------------------------------------------------------------------
# The capture intervals (GRCh38).  Both NKC sub-clusters, each padded by 50 kb:
#
#   chr12:7,670,000-9,220,000   extended  12p13.31  CLEC4C .. KLRG1
#   chr12:9,540,000-10,510,000  core      12p13.2   KLRB1  .. KLRC1
#
# Identical to scripts/srwgs_extract_reads.sh.  They are also the spans the
# decoy contigs in results/srwgs/panel/nkc_refs_decoy.fasta were cut from, so
# changing them here without rebuilding the panel would leave reads with no
# decoy to fall onto.
# --------------------------------------------------------------------------

CONTIG = "chr12"
INTERVALS = [(7_670_000, 9_220_000), (9_540_000, 10_510_000)]

# Public, free to read, and the reference family All of Us aligns against.  The
# M5 in the CRAM header is what actually decides; this is only where the bytes
# come from.
PUBLIC_HG38 = ("gs://gcp-public-data--broad-references/hg38/v0/"
               "Homo_sapiens_assembly38.fasta")
EBI_REF = "http://www.ebi.ac.uk/ena/cram/md5/%s"

FIELDS = ["sample", "status", "elapsed_s", "n_bam",
          "bytes_bam", "bytes_fq", "attempts", "detail"]

OK_STATUS = ("ok",)


# --------------------------------------------------------------------------
# GCS access tokens
# --------------------------------------------------------------------------

class Token:
    """A GCP access token, refreshed before it can expire mid-run.

    Deliberately duplicated from klrc2_cn.py rather than imported: both files
    are meant to be copied into the Workbench on their own.
    """

    CMDS = [
        ["gcloud", "auth", "application-default", "print-access-token"],
        ["gcloud", "auth", "print-access-token"],
    ]

    def __init__(self, ttl: int = 1800, static: str | None = None):
        self.ttl = ttl
        self.static = static or os.environ.get("GCS_OAUTH_TOKEN") or None
        self._val: str | None = None
        self._at = 0.0
        self._lock = threading.Lock()

    def get(self) -> str:
        if self.static:
            return self.static
        with self._lock:
            now = time.monotonic()
            if self._val is not None and now - self._at < self.ttl:
                return self._val
            last = None
            for cmd in self.CMDS:
                try:
                    r = subprocess.run(cmd, capture_output=True, text=True,
                                       timeout=120)
                    if r.returncode == 0 and r.stdout.strip():
                        self._val, self._at = r.stdout.strip(), now
                        return self._val
                    last = r.stderr.strip()
                except (OSError, subprocess.TimeoutExpired) as exc:
                    last = str(exc)
            raise RuntimeError(
                "could not mint a GCP access token. Run `gcloud auth login` (or "
                f"pass --token) and retry. Last error: {last}")


def ref_pattern(cache: str) -> str:
    return os.path.join(cache, "%2s", "%2s", "%s")


def build_env(args, token: Token | None, url: str) -> dict:
    """Environment for one samtools call: credentials plus reference policy."""
    env = dict(os.environ)
    if url.startswith("gs://") and token is not None:
        env["GCS_OAUTH_TOKEN"] = token.get()
    if args.user_project:
        env["GCS_REQUESTER_PAYS_PROJECT"] = args.user_project

    # Unlike read counting, CRAM->BAM decodes SEQ and needs the reference. An
    # unpinned htslib would go looking for it at the EBI registry, which from a
    # locked-down VM hangs rather than fails. Point it at the local M5 cache
    # and nowhere else.
    if args.ref_cache:
        pat = ref_pattern(os.path.abspath(args.ref_cache))
        env["REF_CACHE"] = pat
        env["REF_PATH"] = pat + (":" + EBI_REF if args.allow_remote_ref else "")
    elif not args.allow_remote_ref:
        env["REF_PATH"] = os.devnull
    return env


# --------------------------------------------------------------------------
# Remote object helpers
# --------------------------------------------------------------------------

def gs_to_https(url: str) -> str:
    if url.startswith("gs://"):
        return "https://storage.googleapis.com/" + url[5:]
    return url


def http_get(url: str, token: Token | None, user_project: str,
             byte_range: tuple[int, int] | None = None, timeout: int = 300):
    req = urllib.request.Request(gs_to_https(url))
    unauth = ""
    if token is not None and url.startswith("gs://"):
        try:
            req.add_header("Authorization", f"Bearer {token.get()}")
        except RuntimeError as exc:
            # The default reference bucket is public, and a machine with no
            # gcloud can still read it. Only mention the missing credential if
            # the request actually comes back unauthorised.
            unauth = str(exc)
    if user_project:
        req.add_header("x-goog-user-project", user_project)
    if byte_range:
        req.add_header("Range", f"bytes={byte_range[0]}-{byte_range[1]}")
    try:
        return urllib.request.urlopen(req, timeout=timeout)
    except urllib.error.HTTPError as exc:
        note = (f"\n  no credentials were attached: {unauth}"
                if unauth and exc.code in (401, 403) else "")
        raise RuntimeError(f"{url}: HTTP {exc.code} {exc.reason}{note}") from None
    except urllib.error.URLError as exc:
        raise RuntimeError(f"{url}: {exc.reason}") from None


def read_header(cram: str, crai: str, args, token: Token | None) -> list[str]:
    """@SQ lines of a CRAM header. Needs credentials, no index, no reference."""
    cmd = ["samtools", "view", "-H", "--no-PG"]
    if crai:
        cmd += ["-X", cram, crai]
    else:
        cmd += [cram]
    env = build_env(args, token, cram)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env,
                       timeout=args.timeout)
    if r.returncode != 0:
        raise RuntimeError(f"could not read the CRAM header: "
                           f"{r.stderr.strip().splitlines()[-1] if r.stderr.strip() else 'no error text'}")
    return [l for l in r.stdout.splitlines() if l.startswith("@SQ")]


def sq_tags(sq_lines: list[str]) -> dict[str, dict[str, str]]:
    out = {}
    for line in sq_lines:
        tags = {}
        for field in line.split("\t")[1:]:
            k, _, v = field.partition(":")
            tags[k] = v
        if "SN" in tags:
            out[tags["SN"]] = tags
    return out


def resolve_contig(tags: dict[str, dict[str, str]], want: str) -> str:
    """Pick the header's name for chr12, tolerating an unprefixed build."""
    for name in (want, want.removeprefix("chr"), "chr" + want.removeprefix("chr")):
        if name in tags:
            return name
    raise RuntimeError(
        f"the CRAM header has no contig named {want!r} (nor {want.removeprefix('chr')!r}). "
        f"First few contigs: {list(tags)[:5]}")


def regions_for(contig: str) -> list[str]:
    return [f"{contig}:{a}-{b}" for a, b in INTERVALS]


# --------------------------------------------------------------------------
# The reference cache
# --------------------------------------------------------------------------

def cache_path(cache: str, m5: str) -> str:
    return os.path.join(os.path.abspath(cache), m5[:2], m5[2:4], m5[4:])


def md5_of_file(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def fai_entry(source: str, contig: str, token, user_project) -> tuple[int, int, int, int]:
    """(length, offset, linebases, linewidth) for *contig* from source.fai."""
    if os.path.exists(source):
        text = open(source + ".fai").read()
    else:
        with http_get(source + ".fai", token, user_project) as resp:
            text = resp.read().decode()
    for line in text.splitlines():
        f = line.split("\t")
        if f and f[0] == contig:
            return int(f[1]), int(f[2]), int(f[3]), int(f[4])
    raise RuntimeError(f"{source}.fai has no entry for {contig}")


def fetch_contig(source: str, contig: str, out_path: str, token, user_project,
                 progress=True) -> str:
    """Write the bare sequence of *contig* to out_path; return its MD5.

    A REF_CACHE entry is the sequence with no header, no newlines and no lower
    case — which is exactly what the M5 tag is the checksum of, so writing the
    file and computing the key it must be stored under is one pass.
    """
    length, offset, linebases, linewidth = fai_entry(source, contig, token,
                                                     user_project)
    h = hashlib.md5()
    written = 0
    t0 = time.monotonic()

    if os.path.exists(source):
        # A local FASTA is already free to read; let faidx do the seeking.
        proc = subprocess.Popen(["samtools", "faidx", source, contig],
                                stdout=subprocess.PIPE)
        stream, close = proc.stdout, lambda: proc.wait()
        stream.readline()  # drop the '>' line
    else:
        nlines = -(-length // linebases)
        nbytes = length + nlines * (linewidth - linebases)
        resp = http_get(source, token, user_project,
                        byte_range=(offset, offset + nbytes - 1))
        stream, close = resp, resp.close

    tmp = out_path + f".part{os.getpid()}"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    try:
        with open(tmp, "wb") as fh:
            while written < length:
                buf = stream.read(1 << 23)
                if not buf:
                    break
                buf = buf.translate(None, b"\r\n").upper()
                if written + len(buf) > length:
                    buf = buf[:length - written]   # never run into the next contig
                fh.write(buf)
                h.update(buf)
                written += len(buf)
                if progress:
                    print(f"\r  {written / 1e6:7.1f} / {length / 1e6:.1f} MB",
                          end="", file=sys.stderr, flush=True)
        if progress:
            print(f"   ({time.monotonic() - t0:.0f}s)", file=sys.stderr)
        if written != length:
            raise RuntimeError(f"short read: {written} of {length} bases")
        os.replace(tmp, out_path)
    finally:
        close()
        if os.path.exists(tmp):
            os.unlink(tmp)
    return h.hexdigest()


def ensure_reference(args, token, cram: str, crai: str, build: bool) -> tuple[str, str]:
    """Make sure chr12 is in the cache under the M5 this CRAM asks for.

    Returns (contig name as the header spells it, M5).
    """
    tags = sq_tags(read_header(cram, crai, args, token))
    contig = resolve_contig(tags, args.contig)
    m5 = tags[contig].get("M5", "")
    ur = tags[contig].get("UR", "")

    if not m5:
        raise SystemExit(
            f"the CRAM header carries no M5 for {contig}, so htslib cannot look the "
            "reference up by checksum. Pass --reference <local GRCh38 fasta> instead"
            + (f"; the header points at {ur}" if ur else "") + ".")

    if args.reference:
        return contig, m5

    dest = cache_path(args.ref_cache, m5)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return contig, m5
    if not build:
        raise SystemExit(
            f"no reference for {contig} (M5 {m5}) in {os.path.abspath(args.ref_cache)}.\n"
            f"Build it once with:\n"
            f"    {os.path.basename(sys.argv[0])} ref --from-cram {cram}"
            + (f" \\\n        --crai {crai}" if crai else "")
            + (f"\n(the CRAM header points at {ur})" if ur else ""))

    print(f"fetching {contig} (M5 {m5}) from {args.ref_source}", file=sys.stderr)
    got = fetch_contig(args.ref_source, contig, dest, token, args.user_project)
    if got != m5:
        os.unlink(dest)
        raise SystemExit(
            f"reference mismatch: {contig} from {args.ref_source} has M5 {got}, but the "
            f"CRAM was compressed against M5 {m5}. Decoding against it would corrupt "
            f"every read.\nPoint --ref-source at the reference this cohort actually "
            f"used" + (f" — the header names {ur}" if ur else "") + ".")
    print(f"cached {dest}", file=sys.stderr)
    return contig, m5


# --------------------------------------------------------------------------
# Running samtools
# --------------------------------------------------------------------------

def _fmt(cmd: list[str]) -> str:
    return " ".join(c if not any(ch in c for ch in " *?") else repr(c) for c in cmd)


# htslib reports a 404 as ENOENT. A path that is simply wrong will still be
# wrong on the fourth attempt, and at cohort scale a bad manifest column would
# otherwise burn four timeouts per sample.
PERMANENT = ("no such file or directory", "does not exist", "file not found")


def run_retry(cmd: list[str], args, env_fn, cwd: str) -> tuple[bool, str, int]:
    """Run *cmd*, retrying transient remote failures.

    Object stores reset connections and time out when many workers pull at
    once; on the 1000 Genomes run 6 of 36 samples failed a first pass purely
    from that and all 6 succeeded on retry. A sample dropped for a network
    error is indistinguishable downstream from one that could not be typed.

    The environment is rebuilt per attempt so that a run whose retries span
    longer than the token's life picks up a fresh one.
    """
    detail = ""
    for attempt in range(1, args.retries + 2):
        if attempt > 1:
            time.sleep(min(args.backoff * 2 ** (attempt - 2), 60))
        try:
            env = env_fn()
        except RuntimeError as exc:
            return False, str(exc), attempt
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=args.timeout, env=env, cwd=cwd)
        except subprocess.TimeoutExpired:
            detail = f"timed out after {args.timeout}s"
            continue
        except OSError as exc:
            return False, f"could not run {cmd[0]}: {exc}", attempt
        if r.returncode == 0:
            return True, "", attempt
        lines = [l for l in r.stderr.strip().splitlines() if l.strip()]
        detail = lines[-1] if lines else f"exit {r.returncode}"
        if any(p in detail.lower() for p in PERMANENT):
            return False, detail, attempt
    return False, detail, args.retries + 1


def run_pipeline(stages: list[list[str]], out_path: str, env: dict, cwd: str,
                 timeout: int) -> tuple[bool, str]:
    """Run a pipeline with no shell, failing if any stage fails.

    The shell version relies on pipefail; doing it here means every stage's
    status is checked explicitly and a truncated fastq cannot be mistaken for
    a complete one.
    """
    procs, errs = [], []
    prev = None
    try:
        try:
            with open(out_path, "wb") as out:
                for i, cmd in enumerate(stages):
                    last = i == len(stages) - 1
                    err = tempfile.TemporaryFile()
                    errs.append(err)
                    p = subprocess.Popen(
                        cmd, stdin=prev, stdout=(out if last else subprocess.PIPE),
                        stderr=err, env=env, cwd=cwd)
                    if prev is not None:
                        prev.close()      # only the child holds the read end now
                    prev = p.stdout
                    procs.append(p)
                try:
                    procs[-1].wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    return False, f"fastq pipeline timed out after {timeout}s"
                for p in procs:
                    p.wait()
        except OSError as exc:
            return False, f"could not start the fastq pipeline: {exc}"
        finally:
            if prev is not None:
                prev.close()
            for p in procs:
                if p.poll() is None:
                    p.kill()
                    p.wait()

        for cmd, p, err in zip(stages, procs, errs):
            if p.returncode != 0:
                err.seek(0)
                text = err.read().decode(errors="replace").strip().splitlines()
                return False, f"{cmd[0]}: {text[-1] if text else f'exit {p.returncode}'}"
        return True, ""
    finally:
        for e in errs:
            e.close()


# --------------------------------------------------------------------------
# Destinations
# --------------------------------------------------------------------------

class Dest:
    """Where finished files go: a local directory or a bucket prefix."""

    def __init__(self, path: str, args):
        self.remote = path.startswith("gs://")
        self.path = path.rstrip("/") if self.remote else os.path.abspath(path)
        self.args = args
        self.cli = None
        if self.remote:
            self.cli = ("gcloud" if shutil.which("gcloud")
                        else "gsutil" if shutil.which("gsutil") else None)
            if self.cli is None:
                raise SystemExit("a gs:// destination needs gcloud or gsutil on PATH")
        else:
            os.makedirs(self.path, exist_ok=True)

    # Scratch lives inside a local destination so that finishing a sample is a
    # rename rather than an 80 MB copy across filesystems.
    def staging_parent(self) -> str | None:
        return self.args.tmp_dir or (None if self.remote else self.path)

    def _cmd(self, *rest) -> list[str]:
        if self.cli == "gcloud":
            head = ["gcloud", "storage"]
            tail = ([f"--billing-project={self.args.user_project}"]
                    if self.args.user_project else [])
        else:
            head = ["gsutil"] + (["-u", self.args.user_project]
                                 if self.args.user_project else [])
            tail = []
        return head + list(rest) + tail

    def existing(self) -> set[str]:
        """Samples already staged. One listing, not a stat per sample."""
        if not self.remote:
            names = os.listdir(self.path) if os.path.isdir(self.path) else []
            return {n[:-len(".nkc.fq.gz")] for n in names if n.endswith(".nkc.fq.gz")}
        r = subprocess.run(self._cmd("ls", f"{self.path}/*.nkc.fq.gz"),
                           capture_output=True, text=True)
        if r.returncode != 0 and "matched no objects" not in (r.stderr or "").lower() \
                and "not found" not in (r.stderr or "").lower():
            print(f"WARNING: could not list {self.path}; resume will re-do work:\n"
                  f"  {r.stderr.strip()}", file=sys.stderr)
        out = set()
        for line in r.stdout.splitlines():
            base = line.strip().rsplit("/", 1)[-1]
            if base.endswith(".nkc.fq.gz"):
                out.add(base[:-len(".nkc.fq.gz")])
        return out

    def deliver(self, paths: list[str]) -> tuple[bool, str]:
        if not self.remote:
            for p in paths:
                target = os.path.join(self.path, os.path.basename(p))
                try:
                    os.replace(p, target)
                except OSError:
                    shutil.move(p, target)   # --tmp-dir on another filesystem
            return True, ""
        r = subprocess.run(self._cmd("cp", *paths, self.path + "/"),
                           capture_output=True, text=True)
        if r.returncode != 0:
            lines = r.stderr.strip().splitlines()
            return False, f"upload failed: {lines[-1] if lines else 'unknown'}"
        return True, ""


# --------------------------------------------------------------------------
# Extraction
# --------------------------------------------------------------------------

def slice_cmd(sample: str, cram: str, crai: str, bam: str, regions: list[str],
              args) -> list[str]:
    cmd = ["samtools", "view", "-b", "-@", str(args.threads), "-o", bam]
    if args.reference:
        cmd += ["-T", args.reference]
    cmd += (["-X", cram, crai] if crai else [cram]) + regions
    return cmd


def fastq_stages(sample: str, bam: str, wd: str, args) -> list[list[str]]:
    """The fastq pipeline, kept identical to srwgs_extract_reads.sh.

    `-N` forces the /1 /2 suffixes the caller's read names depend on. The
    defaults matter as much as the flags: -F 0x900 drops secondary and
    supplementary records, and `-s /dev/null` drops singletons — reads whose
    mate fell outside the intervals, about 1% of primary records. Both
    exclusions are part of the validated configuration; changing them here
    would make All of Us calls incomparable to the 1000 Genomes ones.
    """
    return [
        ["samtools", "collate", "-Ou", "-@", str(args.threads),
         "-T", os.path.join(wd, "collate"), bam],
        ["samtools", "fastq", "-N", "-0", os.devnull, "-s", os.devnull, "-"],
        ["gzip"],
    ]


def extract_one(sample: str, cram: str, crai: str, args, token: Token | None,
                regions: list[str], dest: Dest) -> dict:
    """Stage one sample. Never raises; failure is reported in `status`."""
    row = {k: "" for k in FIELDS}
    row["sample"] = sample
    t0 = time.monotonic()

    def env_fn():
        return build_env(args, token, cram)

    # One scratch directory per sample: htslib caches the remote .crai in the
    # working directory, and collate's temp files land here too. Left in a
    # shared directory the indexes alone are ~1.3 MB x the cohort.
    wd = tempfile.mkdtemp(prefix=".nkc_tmp_", dir=dest.staging_parent())
    bam = os.path.join(wd, f"{sample}.nkc.bam")
    fq = os.path.join(wd, f"{sample}.nkc.fq.gz")
    try:
        ok, detail, attempts = run_retry(
            slice_cmd(sample, cram, crai, bam, regions, args), args, env_fn, wd)
        row["attempts"] = attempts
        if not ok:
            row["status"] = "auth_error" if "access token" in detail else "slice_failed"
            row["detail"] = detail
            return row

        # A connection dropped mid-write does not always surface as a non-zero
        # exit, and a silently short slice is a silently under-genotyped
        # sample. quickcheck costs milliseconds.
        qc = subprocess.run(["samtools", "quickcheck", "-v", bam],
                            capture_output=True, text=True)
        if qc.returncode != 0:
            row["status"] = "bam_truncated"
            row["detail"] = (qc.stdout + qc.stderr).strip().replace("\n", "; ")
            return row

        n = subprocess.run(["samtools", "view", "-c", bam],
                           capture_output=True, text=True)
        row["n_bam"] = n.stdout.strip() if n.returncode == 0 else ""
        if row["n_bam"] == "0":
            # Not a crash, but never right for a 2.5 Mb autosomal window: it
            # means the regions missed the header's contig naming.
            row["status"], row["detail"] = "empty_slice", f"no reads in {' '.join(regions)}"
            return row

        if args.index:
            subprocess.run(["samtools", "index", "-@", str(args.threads), bam],
                           capture_output=True, text=True)

        # Everything from here reads the local BAM: no credentials, no
        # reference, nothing that can expire.
        ok, detail = run_pipeline(fastq_stages(sample, bam, wd, args), fq,
                                  dict(os.environ), wd, args.timeout)
        if not ok:
            row["status"], row["detail"] = "fastq_failed", detail
            return row

        row["bytes_bam"] = os.path.getsize(bam)
        row["bytes_fq"] = os.path.getsize(fq)
        if row["bytes_fq"] < 1024:
            row["status"], row["detail"] = "fastq_empty", "fastq is essentially empty"
            return row

        keep = [fq] + ([bam] if args.keep_bam else [])
        if args.index and args.keep_bam and os.path.exists(bam + ".bai"):
            keep.append(bam + ".bai")
        ok, detail = dest.deliver(keep)
        if not ok:
            row["status"], row["detail"] = "deliver_failed", detail
            return row

        row["status"] = "ok"
        return row
    finally:
        row["elapsed_s"] = f"{time.monotonic() - t0:.1f}"
        shutil.rmtree(wd, ignore_errors=True)


# --------------------------------------------------------------------------
# Manifests
# --------------------------------------------------------------------------

SAMPLE_COLS = ("sample", "person_id", "research_id", "sample_id", "id")
CRAM_COLS = ("cram", "cram_url", "cram_uri", "cram_path", "wgs_cram", "path")
CRAI_COLS = ("crai", "crai_url", "crai_uri", "cram_index", "cram_index_uri",
             "cram_index_path", "index")


def sample_from_url(url: str) -> str:
    """Best guess at a sample id: wgs_1000000.cram -> 1000000."""
    base = url.rstrip("/").rsplit("/", 1)[-1]
    for suffix in (".cram", ".bam"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base.removeprefix("wgs_")


def parse_records(spec: str) -> list[tuple[str, str, str]]:
    """A manifest path, or literal `sample,cram[,crai]` records split on ';'."""
    if not os.path.exists(spec):
        if "," not in spec and "://" not in spec:
            raise SystemExit(f"{spec}: no such file, and not a sample,cram[,crai] record")
        out = []
        for rec in spec.split(";"):
            f = [x.strip() for x in rec.split(",") if x.strip()]
            if not f:
                continue
            if "://" in f[0]:                       # cram[,crai], id from the name
                out.append((sample_from_url(f[0]), f[0], f[1] if len(f) > 1 else ""))
            elif len(f) >= 2:
                out.append((f[0], f[1], f[2] if len(f) > 2 else ""))
            else:
                raise SystemExit(f"{rec!r}: expected sample,cram[,crai]")
        return out

    with open(spec, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        delim = "," if head.count(",") > head.count("\t") else "\t"
        low = [c.strip().strip('"').lower() for c in head.rstrip("\n").split(delim)]
        headed = any(c in SAMPLE_COLS + CRAM_COLS + CRAI_COLS for c in low)

        if not headed:
            # Exactly the line a Workbench query hands back: sample,cram,crai
            out = []
            for row in csv.reader(fh, delimiter=delim):
                f = [x.strip() for x in row if x.strip()]
                if not f or f[0].startswith("#"):
                    continue
                if "://" in f[0]:
                    out.append((sample_from_url(f[0]), f[0], f[1] if len(f) > 1 else ""))
                elif len(f) >= 2:
                    out.append((f[0], f[1], f[2] if len(f) > 2 else ""))
            if not out:
                raise SystemExit(f"{spec}: no usable rows, and no header row recognised")
            return out

        rdr = csv.DictReader(fh, delimiter=delim)
        cols = {c.strip().lower(): c for c in (rdr.fieldnames or [])}
        s_col = next((cols[k] for k in SAMPLE_COLS if k in cols), None)
        u_col = next((cols[k] for k in CRAM_COLS if k in cols), None)
        i_col = next((cols[k] for k in CRAI_COLS if k in cols), None)
        if not u_col:
            raise SystemExit(f"{spec}: need a cram column; found {rdr.fieldnames}")
        out = []
        for r in rdr:
            u = (r[u_col] or "").strip()
            if not u:
                continue
            s = (r[s_col] or "").strip() if s_col else ""
            out.append((s or sample_from_url(u), u,
                        (r[i_col] or "").strip() if i_col else ""))
        return out


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------

def check_samtools(need_gcs: bool) -> str:
    if shutil.which("samtools") is None:
        raise SystemExit("samtools is not on PATH")
    v = subprocess.run(["samtools", "--version"], capture_output=True, text=True)
    feat = next((l for l in v.stdout.splitlines()
                 if "GCS=" in l or "libcurl=" in l), "")
    if need_gcs and "GCS=yes" not in feat:
        raise SystemExit(
            "this samtools was built without GCS support, so gs:// URLs will not "
            f"open.\n  features: {feat.strip() or 'unknown'}\n"
            "  install one with GCS=yes (conda-forge/bioconda samtools has it).")
    return feat.strip()


def make_token(args) -> Token | None:
    return None if args.token == "none" else Token(ttl=args.token_ttl,
                                                   static=args.token or None)


def cmd_ref(args) -> int:
    check_samtools(need_gcs=args.from_cram.startswith("gs://") if args.from_cram else False)
    args.reference = ""     # this command builds the cache; -T would skip it
    token = make_token(args)

    if args.from_cram:
        contig, m5 = ensure_reference(args, token, args.from_cram, args.crai,
                                      build=True)
    else:
        if not args.m5:
            raise SystemExit("give --from-cram (preferred) or --m5")
        contig, m5 = args.contig, args.m5
        dest = cache_path(args.ref_cache, m5)
        if os.path.exists(dest):
            print(f"already cached: {dest}", file=sys.stderr)
        else:
            got = fetch_contig(args.ref_source, contig, dest, token,
                               args.user_project)
            if got != m5:
                os.unlink(dest)
                raise SystemExit(f"mismatch: fetched {contig} has M5 {got}, wanted {m5}")

    path = cache_path(args.ref_cache, m5)
    print(f"\n{contig}  M5 {m5}\n  {path}  ({os.path.getsize(path) / 1e6:.1f} MB)")

    if args.also_fasta:
        fa = os.path.join(os.path.abspath(args.ref_cache), f"{contig}.fa")
        with open(path, "rb") as src, open(fa, "wb") as out:
            out.write(f">{contig}\n".encode())
            while True:
                buf = src.read(60 * 10000)
                if not buf:
                    break
                for i in range(0, len(buf), 60):
                    out.write(buf[i:i + 60] + b"\n")
        subprocess.run(["samtools", "faidx", fa], check=False)
        print(f"  FASTA for --reference: {fa}")
    print("\nExtraction can now run with no outbound reference lookups.")
    return 0


def cmd_selftest(args) -> int:
    sample = args.sample or sample_from_url(args.cram)
    if args.dry_run:
        wd = "<scratch>"
        bam = f"{wd}/{sample}.nkc.bam"
        print(_fmt(slice_cmd(sample, args.cram, args.crai, bam,
                             regions_for(args.contig), args)))
        for st in fastq_stages(sample, bam, wd, args):
            print("  | " + _fmt(st))
        return 0

    feat = check_samtools(need_gcs=args.cram.startswith("gs://"))
    print(f"samtools: {subprocess.run(['samtools', '--version'], capture_output=True, text=True).stdout.splitlines()[0]}")
    print(f"htslib features: {feat or 'unknown'}")

    token = make_token(args)
    if args.cram.startswith("gs://") and token is not None:
        try:
            token.get()
            print("access token: ok")
        except RuntimeError as exc:
            print(f"FAIL: {exc}", file=sys.stderr)
            return 1
    if args.user_project:
        print(f"billing project: {args.user_project}")

    try:
        contig, m5 = ensure_reference(args, token, args.cram, args.crai,
                                      build=args.build_ref)
    except SystemExit as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    except RuntimeError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    regions = regions_for(contig)
    print(f"reference: {contig} M5 {m5} "
          f"({args.reference or cache_path(args.ref_cache, m5)})")
    print(f"regions:   {'  '.join(regions)}")

    dest = Dest(args.keep or tempfile.mkdtemp(prefix="nkc_selftest_"), args)

    t0 = time.monotonic()
    row = extract_one(sample, args.cram, args.crai, args, token, regions, dest)
    dt = time.monotonic() - t0

    print(f"\nstatus: {row['status']}   elapsed: {dt:.1f}s   "
          f"attempts: {row['attempts']}")
    if row["detail"]:
        print(f"detail: {row['detail']}")
    if row["status"] != "ok":
        print("\nFAIL — the commands that would run are printed by --dry-run",
              file=sys.stderr)
        return 1

    mb = (int(row["bytes_bam"] or 0) + int(row["bytes_fq"] or 0)) / 1e6
    print(f"records in slice: {row['n_bam']}")
    print(f"bam {int(row['bytes_bam']) / 1e6:.1f} MB   "
          f"fq {int(row['bytes_fq']) / 1e6:.1f} MB   ({mb:.0f} MB per sample)")
    print(f"\nwrote into {dest.path}")
    print(f"projected throughput at --jobs {args.jobs}: "
          f"{args.jobs * 60 / max(dt, 1e-9):.0f} samples/min, "
          f"{args.jobs * 3600 / max(dt, 1e-9) * mb / 1e3:.0f} GB/hour of output")
    return 0


def free_gb(path: str) -> float:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def cmd_extract(args) -> int:
    todo = parse_records(args.manifest)
    remote = any(c.startswith("gs://") for _, c, _ in todo)
    check_samtools(need_gcs=remote)

    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        if not 0 <= i < n:
            raise SystemExit(f"--shard {args.shard}: index must be in [0,{n})")
        todo = [t for k, t in enumerate(todo) if k % n == i]
        print(f"shard {i}/{n}: {len(todo)} samples", file=sys.stderr)

    if args.dry_run:
        s, c, i = todo[0]
        bam = f"<scratch>/{s}.nkc.bam"
        print(_fmt(slice_cmd(s, c, i, bam, regions_for(args.contig), args)))
        for st in fastq_stages(s, bam, "<scratch>", args):
            print("  | " + _fmt(st))
        print(f"\n{len(todo)} samples -> {args.dest}")
        return 0

    dest = Dest(args.dest, args)
    token = make_token(args)

    # One header read up front settles credentials, contig naming and the
    # reference for the whole cohort, instead of surfacing as N identical
    # per-sample failures an hour in. A few candidates, so that one unreadable
    # row cannot stop a cohort that is otherwise fine.
    last = None
    for probe in todo[:3]:
        try:
            contig, m5 = ensure_reference(args, token, probe[1], probe[2],
                                          build=args.build_ref)
            break
        except RuntimeError as exc:
            print(f"  preflight on {probe[0]} failed: {exc}", file=sys.stderr)
            last = exc
    else:
        raise SystemExit(f"preflight failed on the first {len(todo[:3])} samples; "
                         f"last error: {last}")
    regions = regions_for(contig)

    done = dest.existing()
    if done:
        print(f"resuming: {len(done)} samples already staged", file=sys.stderr)
    todo = [t for t in todo if t[0] not in done]
    if not todo:
        print("nothing to do", file=sys.stderr)
        return 0

    # ~80 MB per sample retained, plus a slice in flight per worker. A cohort
    # that runs out of disk halfway is worse than one that refuses to start.
    scratch = dest.staging_parent() or tempfile.gettempdir()
    need = (0 if dest.remote else len(todo) * 0.08) + args.jobs * 0.2
    have = free_gb(scratch)
    print(f"{len(todo)} samples -> {dest.path}   "
          f"(needs ~{need:.1f} GB, {have:.1f} GB free on {scratch})",
          file=sys.stderr)
    if have < need and not args.force:
        raise SystemExit("not enough free space; use a gs:// destination, a bigger "
                         "disk, --tmp-dir elsewhere, or --force")

    log = args.log or os.path.join(
        dest.path if not dest.remote else ".", "nkc_extract_log.tsv")
    fresh = not (os.path.exists(log) and os.path.getsize(log) > 0)
    lock = threading.Lock()
    n_ok = 0

    with open(log, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                           extrasaction="ignore")
        if fresh:
            w.writeheader()
            fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = [ex.submit(extract_one, s, c, i, args, token, regions, dest)
                       for s, c, i in todo]
            for k, fut in enumerate(futures, 1):
                row = fut.result()
                with lock:
                    w.writerow(row)
                    fh.flush()          # a killed job keeps its record
                    if row["status"] in OK_STATUS:
                        n_ok += 1
                    else:
                        print(f"  [{row['status']}] {row['sample']}: {row['detail']}",
                              file=sys.stderr)
                    if k % args.report_every == 0 or k == len(todo):
                        print(f"  {k}/{len(todo)}  ok={n_ok}", flush=True,
                              file=sys.stderr)

    print(f"staged {n_ok}/{len(todo)} samples into {dest.path}", file=sys.stderr)
    print(f"log: {log}", file=sys.stderr)
    return 0 if n_ok == len(todo) else 1


def count_fastq(path: str) -> tuple[int | None, str]:
    """Records in a gzipped fastq, or (None, why not)."""
    gz = subprocess.Popen(["gzip", "-dc", path], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    wc = subprocess.Popen(["wc", "-l"], stdin=gz.stdout, stdout=subprocess.PIPE)
    gz.stdout.close()
    out, _ = wc.communicate()
    err = gz.stderr.read().decode(errors="replace").strip()
    gz.stderr.close()
    gz.wait()
    if gz.returncode != 0 or wc.returncode != 0:
        return None, err.splitlines()[-1] if err else f"gzip exit {gz.returncode}"
    return int(out.split()[0]) // 4, ""


def cmd_verify(args) -> int:
    dest = Dest(args.dest, args)
    if dest.remote:
        raise SystemExit("verify reads the files, so point it at a local directory")
    staged = dest.existing()
    print(f"{len(staged)} fastqs in {dest.path}")

    expected = None
    if args.manifest:
        expected = {s for s, _, _ in parse_records(args.manifest)}
        missing = sorted(expected - staged)
        print(f"manifest: {len(expected)} samples, {len(missing)} missing")
        for s in missing[:20]:
            print(f"  missing {s}")
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more")

    bad, checked, kept = [], 0, []
    for s in sorted(staged):
        fq = os.path.join(dest.path, f"{s}.nkc.fq.gz")
        bam = os.path.join(dest.path, f"{s}.nkc.bam")
        if os.path.getsize(fq) < 1024:
            bad.append((s, "fastq essentially empty"))
            continue
        if os.path.exists(bam):
            qc = subprocess.run(["samtools", "quickcheck", "-v", bam],
                                capture_output=True, text=True)
            if qc.returncode != 0:
                bad.append((s, "bam truncated"))
                continue
        if args.full:
            # Counting reads its way through the whole member, so a truncated
            # fastq — the failure that silently under-genotypes a sample rather
            # than crashing — cannot hide.
            n_fq, err = count_fastq(fq)
            if n_fq is None:
                bad.append((s, f"fastq unreadable: {err}"))
                continue
            if os.path.exists(bam):
                n = subprocess.run(["samtools", "view", "-c", bam],
                                   capture_output=True, text=True)
                if n.returncode == 0 and int(n.stdout.strip()) > 0:
                    kept.append((s, n_fq / int(n.stdout.strip())))
        checked += 1

    print(f"{checked} readable, {len(bad)} bad")
    for s, why in bad[:40]:
        print(f"  {s}: {why}")

    if kept:
        # Reads in the fastq over records in the slice. The shortfall is
        # supplementary records plus mates that fell outside the intervals,
        # ~1% on 1000 Genomes; a sample far off that had a different fetch.
        fracs = sorted(f for _, f in kept)
        mid = fracs[len(fracs) // 2]
        print(f"\nfastq/bam read retention: median {mid:.4f}, "
              f"range {fracs[0]:.4f}-{fracs[-1]:.4f}")
        for s, f in sorted(kept, key=lambda t: t[1])[:5]:
            if abs(f - mid) > 0.02:
                print(f"  outlier {s}: {f:.4f}")
    return 0 if not bad and not (expected and expected - staged) else 1


# --------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def gcs_opts(sp):
        sp.add_argument("--user-project",
                        default=(os.environ.get("GOOGLE_PROJECT")
                                 or os.environ.get("GOOGLE_CLOUD_PROJECT", "")),
                        help="billing project for requester-pays buckets; "
                             "defaults to $GOOGLE_PROJECT [%(default)s]")
        sp.add_argument("--token", default="",
                        help="access token to use verbatim, or 'none' to skip "
                             "token handling (default: $GCS_OAUTH_TOKEN, else "
                             "mint via gcloud)")
        sp.add_argument("--token-ttl", type=int, default=1800,
                        help="seconds before the token is re-minted [%(default)s]")

    def ref_opts(sp):
        sp.add_argument("--ref-cache",
                        default=os.environ.get("NKC_REF_CACHE", "nkc_ref_cache"),
                        help="M5-keyed reference cache directory [%(default)s]")
        sp.add_argument("--ref-source", default=PUBLIC_HG38,
                        help="FASTA the contig is cut from; gs://, https:// or a "
                             "local path [%(default)s]")
        sp.add_argument("--reference", default="",
                        help="use this local FASTA via samtools -T instead of the "
                             "M5 cache")
        sp.add_argument("--allow-remote-ref", action="store_true",
                        help="permit htslib's default reference lookup, which "
                             "reaches out to the EBI CRAM registry")
        sp.add_argument("--contig", default=CONTIG,
                        help="contig holding the NKC [%(default)s]")

    def run_opts(sp):
        sp.add_argument("--jobs", type=int, default=8,
                        help="concurrent samples [%(default)s]")
        sp.add_argument("--threads", type=int, default=2,
                        help="samtools threads per sample [%(default)s]")
        sp.add_argument("--timeout", type=int, default=1800,
                        help="per-command timeout, seconds [%(default)s]")
        sp.add_argument("--retries", type=int, default=3,
                        help="retries per remote command [%(default)s]")
        sp.add_argument("--backoff", type=float, default=10.0,
                        help="initial retry backoff, seconds [%(default)s]")
        sp.add_argument("--tmp-dir", default="",
                        help="parent for per-sample scratch (default: the "
                             "destination when it is local)")
        sp.add_argument("--no-index", dest="index", action="store_false",
                        help="skip samtools index on the slice")
        sp.add_argument("--no-keep-bam", dest="keep_bam", action="store_false",
                        help="deliver only the fastq (the bam is the input to "
                             "depth QC and klrc2_cn.py, so it is kept by default)")
        sp.add_argument("--dry-run", action="store_true",
                        help="print the exact commands and exit")

    r = sub.add_parser("ref", help="build the local reference cache")
    r.add_argument("--from-cram", default="",
                   help="read the contig's M5 from this CRAM's header (preferred)")
    r.add_argument("--crai", default="", help="index for --from-cram")
    r.add_argument("--m5", default="", help="cache this M5 without a CRAM to ask")
    r.add_argument("--also-fasta", action="store_true",
                   help="also write <contig>.fa+.fai for use with --reference")
    ref_opts(r)
    gcs_opts(r)
    r.add_argument("--timeout", type=int, default=600)
    r.set_defaults(func=cmd_ref, reference="")

    s = sub.add_parser("selftest", help="stage one sample and report timing")
    s.add_argument("--cram", required=True)
    s.add_argument("--crai", default="")
    s.add_argument("--sample", default="", help="id for the output files")
    s.add_argument("--keep", default="",
                   help="write into this directory instead of a temp one")
    s.add_argument("--build-ref", action="store_true",
                   help="build the reference cache if it is missing")
    ref_opts(s)
    gcs_opts(s)
    run_opts(s)
    s.set_defaults(func=cmd_selftest)

    e = sub.add_parser("extract", help="stage a cohort")
    e.add_argument("manifest",
                   help="manifest path, or a literal sample,cram[,crai] record")
    e.add_argument("dest", help="output directory or gs:// prefix")
    e.add_argument("--shard", default="",
                   help="process shard i of n, as i/n, for dsub fan-out")
    e.add_argument("--build-ref", action="store_true",
                   help="build the reference cache if it is missing")
    e.add_argument("--log", default="", help="per-sample TSV log")
    e.add_argument("--report-every", type=int, default=25)
    e.add_argument("--force", action="store_true",
                   help="start even if the disk looks too small")
    ref_opts(e)
    gcs_opts(e)
    run_opts(e)
    e.set_defaults(func=cmd_extract)

    v = sub.add_parser("verify", help="check what landed")
    v.add_argument("dest")
    v.add_argument("--manifest", default="")
    v.add_argument("--full", action="store_true",
                   help="also gunzip-test every fastq")
    gcs_opts(v)
    v.set_defaults(func=cmd_verify, tmp_dir="")

    args = p.parse_args()
    try:
        return args.func(args)
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("interrupted; finished samples are already on disk", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
