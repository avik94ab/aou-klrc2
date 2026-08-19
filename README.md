# NKC / KLRC2 tools for All of Us srWGS

Two standalone tools for reading the natural killer complex (chr12p13) out of
All of Us short-read WGS CRAMs, inside the Researcher Workbench. Each is one
Python file, standard library only, so it can be copied into a notebook, a
terminal, or a dsub image without installing anything.

| tool | what it does | output per sample |
|---|---|---|
| [`nkc_extract.py`](nkc_extract.py) | stages the two NKC intervals (2.52 Mb) as BAM + FASTQ for the 27-locus genotyper | ~78 MB |
| [`klrc2_cn.py`](klrc2_cn.py) | measures KLRC2 (NKG2C) copy number from read depth; stages nothing | one row |

Neither downloads a CRAM. Both fetch the index plus the requested slices over
htslib's native GCS backend, so a sample costs tens of megabytes against a CRAM
of tens of gigabytes.

Full documentation is in [nkc_extract.md](nkc_extract.md) and
[klrc2_cn.md](klrc2_cn.md); both include a short list of things to verify inside
the Workbench before committing to a cohort run.

## Getting it into the Workbench

The Workbench has no route back to the machine this was developed on, which is
what this repository is for. In a Workbench terminal:

```bash
git clone https://github.com/avik94ab/aou-klrc2.git
cd aou-klrc2
python nkc_extract.py --help
```

If the VM cannot reach GitHub, `nkc_extract.paste.sh` is the same script encoded
as 205 lines of gzipped base64 that reconstruct the file and check its md5.
Paste the whole thing into a terminal. Regenerate it after editing the script
with `bash make_paste.sh > nkc_extract.paste.sh`.

## Prerequisites

- `samtools` built with GCS support. `samtools --version` must report `GCS=yes`;
  both tools refuse to start otherwise.
- `gcloud`, for minting access tokens. Tokens expire in about an hour, so both
  tools re-mint on a timer rather than reading one at startup.
- Python 3.8 or newer. No third-party packages.
- `gsutil`, only if you write output to a `gs://` destination.

## Quickstart

```bash
CRAM=gs://vwb-aou-datasets-controlled/pooled/wgs/cram/v8_base/wgs_1000000.cram

# build the reference cache: once per VM or image, ~15 s
python nkc_extract.py ref --from-cram $CRAM

# one sample end to end, with timing — run this before anything else
python nkc_extract.py selftest --cram $CRAM --crai $CRAM.crai --sample 1000000 --keep reads/

# a cohort; resumable, and --shard fans out across dsub tasks
python nkc_extract.py extract manifest.csv gs://$WORKSPACE_BUCKET/nkc/reads --jobs 8

# copy number, either from the CRAM or from the slice you just staged
python klrc2_cn.py selftest --cram $CRAM
```

`manifest.csv` is `sample,cram,crai`. Column names from a Workbench query
(`person_id`, `cram_uri`, `cram_index_uri`) are recognised.

All five loci `klrc2_cn.py` measures fall inside the extracted interval, so copy
number can come off the local BAM in 0.1 s instead of 13–16 s of remote queries.

## Data handling

Extracted reads are individual-level data and cannot leave the Workbench. The
genotyper has to come to the data rather than the other way around: only
aggregate results come back out, subject to the usual reporting rules and
small-cell suppression thresholds.

## What is not here

The 27-locus genotyper and its 62 MB HPRC-derived allele panel are not in this
repository. Those go into the workspace bucket through the Workbench's own
storage browser.

## Provenance

Both files are ports of the pipeline used for the published 1000 Genomes run
(3,202 samples). The extraction port reproduces byte-identical FASTQ and BAM
records against the reads that pipeline staged, verified on HG00096 and HG00097.
That parity is the acceptance test rather than a nicety: All of Us calls are
comparable to the 1000 Genomes calls only if read selection is identical, so
nothing in the FASTQ step is configurable.
