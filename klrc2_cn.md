# KLRC2 copy number in All of Us

`klrc2_cn.py` measures KLRC2 copy number straight from remote CRAMs and stages
nothing. The other tool here, `nkc_extract.py`, stages the two NKC intervals as
BAM + FASTQ for the full 27-locus genotyper and is documented separately in
[nkc_extract.md](nkc_extract.md). [README.md](README.md) covers getting either
of them into the Workbench.

`klrc2_cn.py` is a self-contained port of `scripts/srwgs_klrc2_depth.py` +
`scripts/klrc2_population_freq.py` for the All of Us Researcher Workbench.
Standard library only, one file, three subcommands: `selftest`, `measure`, `call`.

The method is unchanged — reads per base at MAPQ ≥ 20 over KLRC2 divided by the
median of three nearby controls (KLRD1, KLRK1, OLR1), then copy number from the
ratio. Only the index and the requested slices are fetched; a CRAM is never
downloaded.

## What changed, and why

### `gs://` instead of `https://`

htslib ships a native GCS backend (`hfile_gcs.so`) driven by two environment
variables, both confirmed present in htslib 1.24:

| variable | purpose |
|---|---|
| `GCS_OAUTH_TOKEN` | bearer token; set from `gcloud auth print-access-token` |
| `GCS_REQUESTER_PAYS_PROJECT` | billing project, if the bucket bills the reader |

Tokens expire in about an hour and a cohort run lasts longer, so `klrc2_cn.py`
re-mints the token on a timer (`--token-ttl`, default 30 min) rather than reading
it once. Each samtools call is a fresh process, so the current value is handed
over at spawn time; there is no long-lived connection holding a stale credential.

### A trimmed target window

The default target is **chr12:10,431,117–10,434,225** rather than the full KLRC2
CDS span. This is the one change that alters numbers, and it was made from
measurement, not preference.

Per-base depth over the full CDS span in deletion homozygotes is not zero — it
sits at about 10% of diploid. That residual is not scattered mismapping: roughly
70% of the span is *exactly* zero and all of the signal is confined to one block
whose left edge falls at chr12:10,434,275–10,434,312 in every homozygote tested
(37 bp of spread across four independent samples from two superpopulations —
the signature of a single ancestral haplotype, consistent with the manuscript's
breakpoint result). Reads there are KLRC2-derived sequence retained on the
deletion chromosome.

Measured on 36 samples of known copy number, in units of `ratio × 2`:

| truth CN | n | trimmed window | full CDS span |
|---|---|---|---|
| 0 | 10 | **0.000 ± 0.001** | 0.190 ± 0.032 |
| 1 | 12 | 0.962 ± 0.107 | 1.144 ± 0.072 |
| 2 | 14 | 1.961 ± 0.123 | 2.124 ± 0.087 |

On the full span the estimator is affine with a non-zero intercept
(`ratio×2 ≈ 0.19 + 0.97 × CN`), which is why the 1000 Genomes run needed a fitted
recalibration before rounding. On the trimmed window the intercept vanishes and
the map is close to the identity, so copy number 0 is a true zero and all 36
genotypes are recovered correctly with **no calibration at all**.

The trade-off is honest: the window is 3,109 bp instead of 4,870, so the CN1 and
CN2 clusters are slightly noisier. The CN0/CN1 gap is nonetheless wider (0.813 vs
0.777). Both windows are always measured, so the comparison can be repeated on
All of Us data before the trimmed window is trusted there.

### Calibration fitted from the data

There is no assembly truth set inside All of Us. The copy-number clusters are far
apart and heavily populated, so the affine map is identified by the data alone:
find the two dominant modes, call them 1 and 2 copies, then refine by alternating
assignment and weighted least squares. On the 1000 Genomes cohort this reproduces
**all 3,197 published calls exactly**, using no labels. Fitting per sequencing
batch also absorbs batch effects for free.

## Running it

```bash
# 1. one sample end to end, with timing — always run this first
python klrc2_cn.py selftest --cram gs://BUCKET/PATH/SAMPLE.cram

# 2. measure a cohort; resumable, and --shard fans out across dsub/Cromwell tasks
python klrc2_cn.py measure manifest.tsv depth.tsv --jobs 16
python klrc2_cn.py measure manifest.tsv depth.$I.tsv --shard $I/$N --jobs 16

# 3. fit the calibration and assign copy number
python klrc2_cn.py call depth.tsv calls.tsv --hist
```

`manifest.tsv` needs a `sample` column and a `cram` (or `cram_url`) column.

Results are appended as each sample lands, so a killed job keeps its work and
re-running skips what is already done. Samples that failed are retried on resume;
samples that succeeded are not.

## Things this handles that are easy to miss

**A failed query is never a zero.** `None` and `0` are distinct throughout, and
the per-sample `status` column records which happened. At the target locus zero is
the *expected* answer for a deletion homozygote, so conflating the two would
manufacture deletions out of network errors.

**Transient remote failures are retried.** In testing, 6 of 36 samples failed on a
first pass at 18 concurrent workers purely from connection resets, and all 6
succeeded on retry. Without retries that is a 17% spurious loss — and it would not
be random with respect to anything you care about.

**Index files are cleaned up.** htslib caches the remote `.crai` in the working
directory. That is what makes five region queries cost one index fetch, but the
indexes are ~1.3 MB each, so a stable working directory accumulates roughly
**335 GB across 250,000 samples**. Each sample gets its own scratch directory,
removed when it finishes.

**No reference is fetched.** `samtools view -c` does not decode SEQ and needs no
reference — verified by running it with `REF_PATH` pointed at a nonexistent path.
But htslib's *default* `REF_PATH` is the EBI CRAM registry, so an unpinned run
inside a locked-down VM can hang on an outbound fetch instead of failing. The
script pins `REF_PATH` to block that; `--allow-remote-ref` restores the default.

**The ambiguity gate is reported per copy-number class.** It fires unevenly — on
1000 Genomes it drops ~10% of CN2 and ~0% of CN0 — so filtering on it biases the
deletion allele frequency. `call` prints the per-class rate, warns when the rates
diverge, and reports the frequency on the *unfiltered* calls. Use the gate for QC,
not for frequency estimation.

## Verify inside the Workbench before a full run

These could not be checked from outside and are worth five minutes each.

1. **`samtools` is built with GCS support.** `samtools --version` must show
   `GCS=yes`. `selftest` checks this and refuses to proceed otherwise.
2. **The CRAM path source.** Locate the srWGS CRAM manifest for the current data
   release rather than assuming a path; bucket layout changes between releases.
3. **Requester-pays.** If reads fail with a 400, set `--user-project`.
4. **Reference build and ALT handling.** All of Us alignments are GRCh38, but
   confirm the contig naming is `chr12` and check whether the pipeline is
   ALT-aware — ALT contigs change which reads clear MAPQ 20 in a paralogous
   region like the KLRC cluster, which is exactly where this method lives.
5. **Rerun the window comparison.** `measure` records both windows. On a few
   hundred All of Us samples, confirm the CN0 cluster really is at zero on the
   trimmed window there too, and refit the calibration before scaling up.
6. **Reporting rules.** Check small-cell suppression thresholds before publishing
   per-population frequencies.

## Scale

Per-sample latency was 13–16 s against EBI over the public internet for all five
regions. In-region GCS reads should be faster. Throughput is network-bound, so it
scales with worker count: `samples/min ≈ workers × 60 / latency`.

For a few hundred thousand samples, shard with `--shard i/n` across dsub tasks.
The data volume is tiny — an index plus a few slices per sample — so cost is
dominated by VM time rather than storage operations, and the whole cohort is a
small compute bill rather than a large one. Estimate it on a real shard before
committing.
