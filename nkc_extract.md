# NKC region extraction in All of Us

`nkc_extract.py` is a self-contained port of `scripts/srwgs_extract_reads.sh`
for the All of Us Researcher Workbench. Standard library only, one file, four
subcommands: `ref`, `selftest`, `extract`, `verify`.

The method is unchanged. All of Us CRAMs are already aligned to GRCh38, so
capture is a ranged fetch over the two NKC intervals rather than anything like
PING's align-everything-to-the-panel step:

| interval | sub-cluster | genes |
|---|---|---|
| `chr12:7,670,000-9,220,000` | extended, 12p13.31 | CLEC4C .. KLRG1 |
| `chr12:9,540,000-10,510,000` | core, 12p13.2 | KLRB1 .. KLRC1 |

Both are the gene span padded by 50 kb, 2.52 Mb in total. About 78 MB moves per
sample against a CRAM of tens of gigabytes, and the CRAM is never downloaded.

## Getting it into the Workbench

The Workbench is a controlled environment: it cannot reach the machine this was
developed on, and there is no shared filesystem. Clone it from a Workbench
terminal:

```bash
git clone https://github.com/avik94ab/aou-klrc2.git
cd aou-klrc2
```

**If the VM cannot reach GitHub, paste it instead** (no laptop step, works in
Verily and Terra). Open a terminal in the Workbench and paste the entire
contents of `nkc_extract.paste.sh`, 205 lines of gzipped base64 that reconstruct
the file and check its md5. Regenerate it whenever the script changes:

```bash
bash make_paste.sh > nkc_extract.paste.sh
```

**Or write it from a notebook cell.** Paste the file contents into a cell under
`%%writefile nkc_extract.py`. Verbose for 800 lines, but it needs no terminal.

**For the genotyping step later, use the bucket.** The panel
(`results/srwgs/panel/`, 62 MB) and `srwgs_genotype_v2.py` are not in this
repository. Those go via `gs://$WORKSPACE_BUCKET/`, uploaded through the
Workbench's own storage browser.

Whichever route, confirm the transfer before trusting it:

```bash
md5sum nkc_extract.py     # must match 7c8c30edc6dbfc06f5486ebbb33bcf68
python nkc_extract.py --help
```

## Running it

```bash
CRAM=gs://vwb-aou-datasets-controlled/pooled/wgs/cram/v8_base/wgs_1000000.cram

# 0. build the reference cache — once per VM or image, ~15 s
python nkc_extract.py ref --from-cram $CRAM

# 1. one sample end to end, with timing — always run this first
python nkc_extract.py selftest --cram $CRAM --crai $CRAM.crai \
    --sample 1000000 --keep reads/

# 2. the sample you have in hand, straight from its manifest line
python nkc_extract.py extract "1000000,$CRAM,$CRAM.crai" reads/

# 3. a cohort; resumable, and --shard fans out across dsub tasks
python nkc_extract.py extract manifest.csv gs://$WORKSPACE_BUCKET/nkc/reads --jobs 8
python nkc_extract.py extract manifest.csv reads/ --shard $I/$N --jobs 8

# 4. confirm what landed is complete and readable
python nkc_extract.py verify reads/ --manifest manifest.csv --full
```

`manifest.csv` is `sample,cram,crai` — with or without a header row, comma or
tab separated. Column names from a Workbench query (`person_id`, `cram_uri`,
`cram_index_uri`, and others) are recognised. **The sample id is taken verbatim
from the first column**, so `000000,gs://.../wgs_1000000.cram,...` produces
`000000.nkc.*`; the file name is only consulted when there is no id column.

Per sample this writes:

| file | |
|---|---|
| `{sample}.nkc.bam` | the raw slice, ~36 MB |
| `{sample}.nkc.bam.bai` | index; `--no-index` to skip |
| `{sample}.nkc.fq.gz` | reads, names suffixed `/1` `/2`, ~42 MB |

## What changed, and why

### A reference is now mandatory

This is the one difference that will stop a run, and it does not apply to
`klrc2_cn.py`. That script pins `REF_PATH` at `/dev/null` because
`samtools view -c` never decodes SEQ, so it needs no reference at all. Writing
a BAM *does* decode SEQ, and CRAM decompression needs the exact sequence the
CRAM was compressed against — not a compatible build, the same bases.

Letting htslib find it on its own is the bad option twice over: its default
`REF_PATH` is the EBI CRAM registry, so an unpinned run inside a locked-down VM
hangs on an outbound fetch instead of failing, and if it ever did resolve
something, a near-miss reference corrupts every read rather than erroring.

`nkc_extract.py ref` closes both holes. It reads the M5 checksum of chr12 out
of the CRAM's own header, byte-range-fetches that single contig out of a public
GRCh38 FASTA — 133 MB moved, not the 3.1 GB file — and stores it under the M5
htslib will look for. Extraction then runs with `REF_PATH` pointed at that
cache and nowhere else. A reference that does not match fails at cache-build
time, loudly, quoting the reference the CRAM header itself names:

```
reference mismatch: chr12 from gs://gcp-public-data--broad-references/... has
M5 <a>, but the CRAM was compressed against M5 <b>. Decoding against it would
corrupt every read.
Point --ref-source at the reference this cohort actually used — the header
names /path/from/the/UR/tag
```

On the 1000 Genomes CRAMs the check passes: chr12 of the public Broad hg38 is
`96e414eace405d8c27a6d35ba19df56f`, byte-identical to the sequence those CRAMs
were compressed against. All of Us is also GRCh38, so it is expected to pass
there too — but it is verified rather than assumed, and step 0 above is where
you find out.

### `gs://` instead of `https://`

Same mechanism as `klrc2_cn.py`: htslib's native GCS backend, driven by
`GCS_OAUTH_TOKEN` and, when the bucket bills the reader,
`GCS_REQUESTER_PAYS_PROJECT`. Tokens expire in about an hour, so the token is
re-minted on a timer and rebuilt for each retry rather than read once at
startup. `--user-project` defaults to `$GOOGLE_PROJECT`.

### Destinations can be buckets

A Workbench VM or a dsub task has a small disk and this produces 78 MB per
sample. A `gs://` destination stages each sample locally, uploads it, and drops
the local copy; resume is answered from a single bucket listing rather than a
stat per sample. A local destination keeps the shell script's behaviour, with
scratch inside the destination so that finishing a sample is a rename rather
than an 80 MB copy.

## Parity, measured

The port is not merely equivalent, it is identical. Both samples below were
re-extracted with `nkc_extract.py` and compared against the files the Wynton
pipeline staged for the published cohort:

| | HG00096 | HG00097 |
|---|---|---|
| fastq, decompressed md5 | identical | identical |
| bam alignment records md5 | identical | identical |

Only the `@PG` provenance lines differ, since they record the command line.
Reproduce it with:

```bash
zcat old/HG00096.nkc.fq.gz | md5sum;  zcat new/HG00096.nkc.fq.gz | md5sum
samtools view old/HG00096.nkc.bam | md5sum;  samtools view new/HG00096.nkc.bam | md5sum
```

That matters because All of Us calls are only comparable to the 1000 Genomes
calls if the reads reaching the caller were selected the same way. Anything in
the fastq step is therefore fixed, not configurable.

## The 1% that is dropped on purpose

`samtools fastq -N -0 /dev/null -s /dev/null` keeps fewer records than the
slice contains, in two ways that are easy to misread. For HG00096, of 653,604
records in the slice:

| | records | |
|---|---|---|
| supplementary | 1,771 | dropped by the `-F 0x900` default |
| singletons | 6,257 | dropped by `-s /dev/null` — mate outside the intervals |
| written to the fastq | 645,576 | 98.77% |

The comment at `scripts/srwgs_extract_reads.sh:63` says these singletons "all go
in one file". They do not; `-s` discards them. The behaviour is fine — it is
what every published call was made from — but the reason to know it is that
98.77% is the number to expect, and `verify --full` reports exactly this
retention so a sample fetched some other way stands out.

## Things this handles that are easy to miss

**A truncated slice is not always a failed command.** A connection dropped
mid-write can leave a short BAM behind a zero exit status, and a short slice is
a silently under-genotyped sample rather than a crash. Every slice is run
through `samtools quickcheck` before it is accepted, and every stage of the
fastq pipeline has its exit status checked individually rather than relying on
`pipefail`.

**An empty slice is caught.** Zero reads over 2.5 Mb of autosome is never
right; it means the regions missed the header's contig naming. That is reported
as `empty_slice`, not as success. The contig name is resolved from the header
anyway, so `12` works as well as `chr12`.

**A wrong path is not retried.** htslib reports a 404 as `No such file or
directory`. Transient resets are retried with backoff — on the 1000 Genomes run
6 of 36 samples failed a first pass purely from concurrency and all 6 succeeded
on retry — but a path that is simply wrong fails immediately, instead of
burning four timeouts per row on a mistyped manifest column.

**One bad row does not stop a cohort.** Preflight tries up to three manifest
rows before concluding that credentials or the reference are the problem.

**Index files are cleaned up.** htslib caches the remote `.crai` in the working
directory. Each sample gets its own scratch directory, removed when it
finishes; left in a shared directory the indexes alone are ~1.3 MB per sample.

**The disk is checked before the run, not during it.** A cohort that runs out
of space halfway is worse than one that refuses to start.

**Extracted reads are individual-level data.** They cannot leave the Workbench.
The genotyper has to come to the data: `scripts/srwgs_genotype_v2.py` and
`results/srwgs/panel/` need to be uploaded into the workspace, and only
aggregate results come back out, subject to the usual reporting rules.

## Verify inside the Workbench before a full run

These could not be checked from outside, and are worth five minutes each.

1. **`samtools` is built with GCS support.** `samtools --version` must show
   `GCS=yes`; `selftest` refuses to proceed otherwise. The environment used on
   Wynton (htslib 1.24) has it.
2. **The reference matches.** `ref --from-cram` is the whole test. If it
   reports a mismatch, the `UR:` tag it prints names the reference All of Us
   actually used; point `--ref-source` at that.
3. **Requester-pays.** If reads fail with a 400, set `--user-project`.
4. **The CRAM path source.** Confirm the manifest for the current data release
   rather than assuming a bucket layout; it changes between releases.
5. **Timing.** Run `selftest` and read the projected throughput off it before
   committing to a shard count.

## Scale

Measured from outside GCP, against 1000 Genomes CRAMs over HTTPS: **37–40 s per
sample** at `--threads 4`, 78 MB of output. In-region GCS reads should be
faster. Throughput is network-bound, so it scales with `--jobs`:
`samples/min ≈ jobs × 60 / latency`, about 13 samples/min at `--jobs 8`.

Storage is the constraint that bites, not compute: 78 MB per sample is 7.8 TB
per 100,000 samples. Use `--no-keep-bam` if the depth QC is not wanted (that
halves it), and write to a bucket rather than a VM disk. Price a real shard
before committing.

For fan-out, `--shard i/n` splits the manifest deterministically:

```bash
dsub ... --command 'python nkc_extract.py extract ${MANIFEST} ${OUT} \
    --shard ${I}/${N} --jobs 8 --build-ref' \
    --env I=0 --env N=200 ...
```

`--build-ref` lets each task construct its own reference cache, which costs 15 s
and 133 MB of local disk per task and avoids shipping a reference into the
image.

## Reusing the slice

All five loci `klrc2_cn.py` measures — KLRC2, KLRD1, KLRK1, OLR1 — fall inside
the core interval, so copy number can be measured from the local BAM instead of
re-querying the CRAM. Measured on HG00096: **0.1 s against the slice versus
13–16 s remote**, ratio 0.9938, i.e. 2 copies, matching the remote answer.

```bash
printf 'sample\tcram\n1000000\treads/1000000.nkc.bam\n' > local.tsv
python klrc2_cn.py measure local.tsv depth.tsv --token none --jobs 8
```

This needs the `.bai`, which is why it is written by default.
