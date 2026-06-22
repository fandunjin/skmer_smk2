# skmer-smk2

`skmer-smk2` is a packaged Snakemake workflow for building phylogenetic trees
from paired-end FASTQ reads. It provides one command for read cleaning, optional
reference-based plastid filtering, base-aware depth normalization, and tree
inference with Skmer, WASTER, and Mash.

## Recommended Running Modes

Most users should start with one of these modes. The detailed sections below
explain the same commands more fully.

### Mode 1: Standard Full Run

Use this when the dataset is moderate in size and you want Skmer, WASTER, and
Mash from one command.

With plastid/reference removal:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2 \
  --fastp-threads 4 \
  --repair-threads 2 \
  --bbmerge-threads 2 \
  --total-mem-mb 180000 \
  --bowtie2-mem-mb 6000 \
  --printshellcmds
```

Without plastid/reference removal:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 --printshellcmds
```

### Mode 2: Recommended Large-Dataset HPC Run

Use this when there are many samples or when WASTER is slow or memory-heavy.
The key idea is to run shared preprocessing first, then run Skmer, Mash, and
WASTER as separate scheduler jobs in the same `RUN_DIR`.

Export editable scripts:

```bash
skmer-smk2 init -o skmer_smk2_templates
```

Edit `FASTQ_DIR`, optional `REF_FASTA`, `RUN_DIR`, scheduler headers, and
environment activation in these files:

```text
stage_preprocess.sh
stage_skmer.sh
stage_mash.sh
stage_waster.sh
```

Then submit:

```bash
jsub < skmer_smk2_templates/stage_preprocess.sh
jsub < skmer_smk2_templates/stage_skmer.sh
jsub < skmer_smk2_templates/stage_mash.sh
jsub < skmer_smk2_templates/stage_waster.sh
```

Use the same `RUN_DIR` in all four scripts so Snakemake can reuse existing
outputs. WASTER is usually the least parallel final stage, so the WASTER script
requests few cores but high memory.

### Mode 3: Check And Repair Inputs First

Use this when FASTQ files may be truncated, gzip-damaged, or have R1/R2 pairing
problems. This is the safest mode before the main analysis.

```bash
skmer-smk2 input-check -i RAW_FASTQ_DIR -o input_qc -j 12
skmer-smk2 input-repair -i RAW_FASTQ_DIR -o input_qc --samples "sampleA sampleB" -j 4
skmer-smk2 run -i input_qc/input_for_skmer -ref REF_FASTA -s 75 -j 48 \
  --bowtie2-threads 2 --fastp-threads 4 --repair-threads 2 --bbmerge-threads 2
```

If no samples need repair, run the final command directly on the original
FASTQ directory or on `input_qc/input_for_skmer/` after confirming the report.

### Mode 4: Resume Or Run One Branch Only

Snakemake reuses existing outputs. After an interrupted run, or when only one
result branch is needed, use the branch flags:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -prep
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 16 -skmer
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -mash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 1 -waster \
  --waster-threads 1 --waster-mem-mb 220000 --total-mem-mb 240000
```

`-prep` runs only shared preprocessing and base-aware normalization. `-skmer`,
`-mash`, and `-waster` run only the selected analysis branch. If none of these
flags is supplied, all three analysis branches run.

### Quick Decision Guide

| Situation | Recommended mode |
| --- | --- |
| Clean moderate dataset, want all outputs | Mode 1 |
| Many samples on HPC, want better resource use | Mode 2 |
| FASTQ damage or uncertain input quality | Mode 3 |
| WASTER is using one core while a large job is allocated | Mode 2 or Mode 4 `-waster` |
| Continue after a failed or killed run | Mode 4 |

The same `skmer-smk2` commands can be used on a local workstation or inside HPC
scheduler scripts. HPC does not use a different workflow; the scheduler only
submits these commands to compute nodes.

## What This Workflow Provides

This package implements the requested runtime-controlled workflow:

| Requirement | Implementation |
| --- | --- |
| `-i` input FASTQ directory | `skmer-smk2 run -i FASTQ_DIR` |
| Optional plastid/reference removal | `-ref REF_FASTA` enables `bowtie2-build`, `bowtie2`, and pair repair; omitting `-ref` skips those rules |
| Statistics after optional `-ref` filtering | `results/stats/post_filter_summary.sorted.tsv` records reads, total bases, and average length after cleaning/filtering/merging |
| Base-aware head instead of fixed read count | `-s 75` selects the base count at the 75% position after sorting samples by total bases, then writes FASTQ records until that base cutoff is reached |
| Threshold selection report | `results/stats/head_cutoff_candidates.tsv` compares candidate cutoffs such as 50, 75, and 90 |
| Head statistics | `results/stats/head_summary.sorted.tsv` records final reads, total bases, and average length for each normalized sample |
| Skmer trees | direct tree, bootstrap consensus tree, and merged tree |
| WASTER trees | WASTER topology tree and branch-length tree |
| Mash outputs | direct tree, bootstrap consensus tree, merged tree, and distance heatmap |

## Workflow Overview

For each paired-end sample, the workflow runs these steps:

1. Detect complete R1/R2 FASTQ pairs from the input directory.
2. Clean reads with `fastp`.
3. If `-ref` is supplied, build a Bowtie2 index and remove reads mapping to the
   reference genome.
4. If `-ref` is supplied, repair synchronized pairs with BBMap `repair.sh`.
5. Merge overlapping reads with BBMap `bbmerge.sh`.
6. Concatenate merged and unmerged reads into one FASTQ per sample.
7. Count reads, total bases, and average read length.
8. Sort samples by total bases from high to low and select the `-s` percentile
   base cutoff.
9. Write normalized FASTQ files by accumulating complete FASTQ records until the
   selected base cutoff is reached.
10. Build Skmer, WASTER, and Mash outputs from the normalized FASTQ files.

## Installation

Install directly from GitHub:

```bash
python -m pip install git+https://github.com/fandunjin/skmer_smk2.git
```

If `git clone` over HTTPS is unstable on your system, install from the GitHub ZIP
archive instead:

```bash
python -m pip install --no-cache-dir --force-reinstall \
  https://github.com/fandunjin/skmer_smk2/archive/refs/heads/main.zip
```

Install from a local clone:

```bash
git clone https://github.com/fandunjin/skmer_smk2.git
cd skmer_smk2
python -m pip install .
```

Check that the command is available:

```bash
skmer-smk2 --version
skmer-smk2 -h
skmer-smk2 run -h
```

## External Software

The Python package installs the command wrapper and the bundled Snakemake
workflow. The bioinformatics programs used by the workflow must also be
available in the active environment.

Required tools:

```text
snakemake
fastp
bowtie2
bowtie2-build
repair.sh
bbmerge.sh
skmer
fastme
raxmlHPC
waster
waster_branchlength
mash
jellyfish
seqtk
seqkit
gzip
```

Check the current environment:

```bash
skmer-smk2 doctor
```

`doctor` is advisory by default. It prints tools found in `PATH`, compatible
tool aliases, and installation suggestions. Use strict mode when missing tools
should make the command fail:

```bash
skmer-smk2 doctor --strict
```

Install conda-available missing tools into the active environment:

```bash
skmer-smk2 doctor --install
```

The automatic installer uses `mamba` if available, otherwise `conda`. It only
installs packages that are available through conda channels. Tools such as
`skmer`, `waster`, and `waster_branchlength` may need manual installation
depending on your environment.

A typical conda/mamba command for common dependencies is:

```bash
mamba install -c conda-forge -c bioconda \
  snakemake fastp bowtie2 bbmap fastme raxml mash seqkit gzip
```

## Input FASTQ Files

`-i` must point directly to the directory containing FASTQ files, not to the
repository root.

Supported paired-end naming styles:

```text
SampleA_1.fq.gz      SampleA_2.fq.gz
SampleB_R1.fq.gz     SampleB_R2.fq.gz
SampleC.R1.fq.gz     SampleC.R2.fq.gz
SampleD-R1.fq.gz     SampleD-R2.fq.gz
```

The sample name is the shared prefix before the mate suffix. For example,
`SampleC.R1.fq.gz` and `SampleC.R2.fq.gz` are detected as sample `SampleC`.

Supported FASTQ extensions:

```text
.fq.gz
.fastq.gz
.fq
.fastq
```

## Input Data Check And Repair

Use two explicit steps before the main analysis: first check every input FASTQ
and review the report, then repair only the sample names you choose.

Step 1: check input data.

```bash
skmer-smk2 input-check -i /path/to/raw_fastq -o input_qc -j 12
```

`-j/--jobs` controls how many FASTQ files are checked at the same time. The
default is `1` to be gentle on shared filesystems. On HPC nodes, `-j 8` to
`-j 16` is usually a practical starting point.

Main check outputs:

```text
input_qc/input_file_report.tsv
input_qc/input_sample_report.tsv
input_qc/logs/*.gzip.log
```

`input_file_report.tsv` records each file's sample name, mate, path, file size,
gzip status, seqkit status, reads, bases, and average length.

`input_sample_report.tsv` is the main list to inspect. It contains:

```text
sample
r1_file
r2_file
r1_size_bytes
r2_size_bytes
r1_gzip_status
r2_gzip_status
r1_reads
r2_reads
r1_bases
r2_bases
r1_avg_len
r2_avg_len
pair_status
suggested_action
note
```

Suggested actions:

| Action | Meaning |
| --- | --- |
| `USE_AS_IS` | Original R1/R2 passed checks and can be used directly |
| `REPAIR_CANDIDATE` | The sample should be considered for `seqkit sana` plus `repair.sh` |
| `MISSING_PAIR` | R1 or R2 is missing |
| `CHECK_MANUALLY` | Information is unusual and should be inspected manually |

Step 2: repair selected samples by name.

```bash
skmer-smk2 input-repair -i /path/to/raw_fastq -o input_qc \
  --samples "sampleA sampleB" -j 4 --repair-threads 2
```

Comma-separated names are also accepted:

```bash
skmer-smk2 input-repair -i /path/to/raw_fastq -o input_qc \
  --samples sampleA,sampleB -j 4 --repair-threads 2
```

For `input-repair`, `-j/--jobs` controls how many selected samples are repaired
or linked at the same time. `--repair-threads` controls how many threads each
`repair.sh` process uses.

Repair outputs:

```text
input_qc/post_repair_file_report.tsv
input_qc/post_repair_sample_report.tsv
input_qc/repaired/
input_qc/input_for_skmer/
```

Only samples listed in `--samples` are repaired. Normal samples not listed in
`--samples` are linked into `input_for_skmer/`. Abnormal samples not listed in
`--samples` are not included, which avoids accidentally analyzing known bad
inputs.

After reviewing `post_repair_sample_report.tsv`, use
`input_qc/input_for_skmer/` as the FASTQ input directory for the main workflow.

Optional helper scripts:

```bash
skmer-smk2 repair-fastq --workdir input_qc --copy-only
```

This writes:

```text
input_qc/01_input_check.sh
input_qc/02_input_repair.sh
```

Edit only the scheduler header and environment activation lines if using them
on an HPC system.

## Running The Workflow

Run the main analysis on the checked/repaired input directory.

With reference filtering:

```bash
skmer-smk2 run -i input_qc/input_for_skmer -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2 \
  --fastp-threads 4 \
  --repair-threads 2 \
  --bbmerge-threads 2
```

Without reference filtering:

```bash
skmer-smk2 run -i input_qc/input_for_skmer -s 75 -j 48
```

You can still run directly on a raw FASTQ directory if you intentionally skip
the input-check step:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2 \
  --fastp-threads 4 \
  --repair-threads 2 \
  --bbmerge-threads 2 \
  --total-mem-mb 180000 \
  --bowtie2-mem-mb 6000
```

Run only selected analysis branches:

Shared preprocessing only:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -prep
```

Mash only:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -mash
```

Skmer and Mash, without WASTER:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -skmer -mash
```

WASTER only:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 1 -waster \
  --waster-threads 1 --waster-mem-mb 220000 --total-mem-mb 240000
```

If none of `-prep`, `-skmer`, `-waster`, or `-mash` is supplied, all three
analysis branches run after preprocessing. `-prep` runs only the shared
preprocessing and base-aware head outputs.

Preview the workflow without running jobs:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 1 --dry-run
```

Print shell commands while running:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 --printshellcmds
```

Pass extra Snakemake arguments after `--`:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 -- --keep-going
```

## Main Runtime Options

| Option | Meaning |
| --- | --- |
| `-i`, `--input` | Directory containing paired FASTQ files |
| `-ref`, `--ref` | Optional reference FASTA used for Bowtie2 filtering |
| `-s`, `--sample-percentile` | Percentile position used to choose the base-count cutoff; default `75` |
| `--candidate-percentiles` | Percentiles included in the cutoff selection report; default `50,60,70,75,80,90,95` |
| `-prep` | Run only shared preprocessing, filtering, merging/fallback, statistics, and `nDNAOK` outputs |
| `-skmer` | Run the Skmer branch |
| `-waster` | Run the WASTER branch |
| `-mash` | Run the Mash branch |
| `-j`, `--jobs` | Snakemake cores/jobs |
| `-b`, `--bootstraps` | Bootstrap replicate count for Skmer and Mash; default `100` |
| `--exclude-samples` | Comma- or whitespace-separated sample names to skip |
| `--bbmerge-timeout` | Seconds to wait for each BBMerge job before fallback; default `14400`, use `0` to disable |
| `--skmer-sketch-size` | Skmer `-s` sketch size for `reference` and `subsample`; default `100000` |
| `--skmer-threads` | Threads used inside Skmer reference/subsample steps; default `16` |
| `--skmer-mem-mb` | Snakemake memory resource for Skmer reference/subsample jobs; default `120000` |
| `--fastp-threads` | Threads per sample for `fastp`; default `4` |
| `--bowtie2-threads` | Threads per sample for Bowtie2 filtering; default `2` |
| `--repair-threads` | Threads per sample for `repair.sh`; default `2` |
| `--bbmerge-threads` | Threads per sample for `bbmerge.sh`; default `2` |
| `--waster-threads` | Threads used by WASTER and `waster_branchlength`; default `4` |
| `--fastp-mem-mb` | Snakemake memory resource per `fastp` job; default `2000` |
| `--bowtie2-mem-mb` | Snakemake memory resource per Bowtie2 filtering job; default `4000` |
| `--repair-mem-mb` | Snakemake memory resource per `repair.sh` job; default `4000` |
| `--bbmerge-mem-mb` | Snakemake memory resource per `bbmerge.sh` job; default `4000` |
| `--waster-mem-mb` | Snakemake memory resource for WASTER and `waster_branchlength`; default `120000` |
| `--total-mem-mb` | Optional total Snakemake `mem_mb` scheduling limit |
| `--workdir` | Directory where `results/` and workflow cache are written; default current directory |
| `--scheduler` | Snakemake scheduler; default `greedy`, use `ilp` only if CBC/ILP solver is installed |
| `--latency-wait` | Snakemake latency wait seconds; default `120` |
| `--dry-run` | Build and print the DAG without running jobs |
| `--printshellcmds` | Print shell commands from Snakemake |

## Per-Sample Parallelism And Memory Scheduling

`-j/--jobs` is the total Snakemake core budget. The per-step thread options
control how many cores each sample job consumes. This means:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2
```

can schedule roughly 24 Bowtie2 sample jobs at the same time, instead of giving
one sample all 48 cores. This is usually faster for many-sample datasets because
Bowtie2 plastid filtering often does not scale efficiently to very high thread
counts for a single sample.

Recommended high-parallel reference-filtering command:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2 \
  --fastp-threads 4 \
  --repair-threads 2 \
  --bbmerge-threads 2 \
  --total-mem-mb 180000 \
  --bowtie2-mem-mb 6000 \
  --printshellcmds
```

If memory pressure is high, lower the total scheduling limit or increase the
per-job memory estimate so Snakemake runs fewer samples at once:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2 \
  --total-mem-mb 96000 \
  --bowtie2-mem-mb 8000
```

`--total-mem-mb` is a Snakemake scheduling resource. It limits how many jobs are
started together, but it is not a hard memory cap enforced by the operating
system or scheduler. Real memory allocation still depends on the HPC job request
and the programs themselves.

## Staged HPC Runs

For large datasets, it is often more efficient to split the workflow into
separate scheduler jobs while keeping the same `RUN_DIR`. The first job runs
shared preprocessing with many sample-level jobs in parallel. Later jobs reuse
those outputs and run only Skmer, Mash, or WASTER.

Export example scripts:

```bash
skmer-smk2 init -o skmer_smk2_templates
```

The staged templates are:

```text
stage_preprocess.sh
stage_skmer.sh
stage_mash.sh
stage_waster.sh
```

Edit `FASTQ_DIR`, optional `REF_FASTA`, `RUN_DIR`, scheduler headers, and
environment activation in each script. Use the same `RUN_DIR` for every stage,
otherwise Snakemake cannot reuse previous outputs.

Recommended order:

```bash
jsub < stage_preprocess.sh
jsub < stage_skmer.sh
jsub < stage_mash.sh
jsub < stage_waster.sh
```

Skmer and Mash can be submitted after preprocessing finishes. WASTER can also be
submitted after preprocessing finishes, but it is commonly the least parallel
final stage and should usually be requested as a few-core, high-memory job.

### Which Steps Use Few Cores?

These steps can use multiple cores per sample and also run many samples in
parallel through `-j`: `fastp`, optional Bowtie2 filtering, `repair.sh`,
`bbmerge.sh`, and Skmer `reference/subsample`.

These steps are mostly single-process but can still use `-j` effectively when
there are many samples or bootstrap replicates: Mash sketching, Mash bootstrap
sketching, FastME tree building for bootstrap replicates, and small Python
conversion/report steps.

These steps are usually one DAG job at a time near the end of a branch:
`mash paste`, `mash dist`, RAxML consensus, tree merging, WASTER, and
`waster_branchlength`. When only WASTER remains, `-j 48` will not make it use
48 cores. In that case, run `-waster` separately with low scheduler cores and
enough memory, for example:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 \
  -j 1 -waster \
  --waster-threads 1 \
  --waster-mem-mb 220000 \
  --total-mem-mb 240000 \
  --printshellcmds
```

## How `-ref` Works

When `-ref REF_FASTA` is provided:

1. `bowtie2-build` builds an index from the reference FASTA.
2. `bowtie2` aligns cleaned read pairs to the reference.
3. Concordantly mapped pairs are removed.
4. Unmapped pairs are repaired with `repair.sh`.
5. The repaired non-reference reads continue to merging, statistics, and tree
   inference.

When `-ref` is omitted, the Bowtie2 and pair-repair branch is skipped. Cleaned
R1/R2 reads go directly to `bbmerge.sh`.

This makes the reference removal optional at runtime without requiring a
different Snakefile.

## How `-s` Base-Aware Normalization Works

The workflow does not use a fixed read count such as 25,000,000 reads. Instead,
it normalizes by total bases:

1. After cleaning and optional reference removal, each sample is converted into
   `results/<sample>/sample.fq`.
2. The workflow counts reads, total bases, and average length for every sample.
3. Samples are sorted by `total_bases` from high to low.
4. `-s 75` means: choose the sample located at the 75% position in this sorted
   list and use that sample's `total_bases` as the cutoff.
5. Each sample is then written to `results/<sample>/nDNAOK/<sample>.fq` by
   accumulating complete FASTQ records until the cutoff base count is reached.

This is length-aware. A 100 bp sample and a 150 bp sample are compared by total
bases, not by read count alone.

## Choosing A Good `-s` Value

The default `-s 75` is a practical middle choice for many datasets, but some
projects may benefit from a higher or lower cutoff. The workflow writes a
candidate report before final head normalization:

```text
results/stats/head_cutoff_candidates.tsv
```

By default, this report compares:

```text
50,60,70,75,80,90,95
```

Important columns:

| Column | Meaning |
| --- | --- |
| `percentile` | Candidate `-s` value |
| `selected` | `yes` for the value actually used in the current run |
| `position` | Sample position after sorting by total bases from high to low |
| `sample` | Sample at that position |
| `cutoff_bases` | Base count that would be used as the shared head cutoff |
| `samples_truncated` | Number of samples larger than the cutoff and therefore trimmed down |
| `samples_below_cutoff` | Number of samples smaller than the cutoff and therefore kept completely |
| `estimated_retained_percent` | Estimated percentage of total bases retained across all samples |
| `smallest_sample_percent_of_cutoff` | Whether the smallest sample is far below the proposed cutoff |

Use the report like this:

1. Check `post_filter_summary.sorted.tsv` to see the full depth distribution.
2. Check `head_cutoff_candidates.tsv` to compare 50, 75, 90, or custom values.
3. Lower values such as 50 choose a higher cutoff, so they retain more total
   bases from deep samples but may leave low-depth samples below the target.
4. Higher values such as 90 choose a lower cutoff, so they trim more data from
   deep samples but make the normalized dataset more balanced.
5. `samples_below_cutoff` and `smallest_sample_percent_of_cutoff` help decide
   whether the chosen cutoff is too high for the smallest samples.

Custom candidate values can be requested without changing the actual `-s`
selection:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 \
  --candidate-percentiles 40,50,60,75,90
```

## Example Data

The repository includes small demonstration files:

```text
raw_data/raw_data/
raw_data/ref.fna
```

Run a dry run with the bundled demo data:

```bash
skmer-smk2 run -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 1 -b 2 --dry-run
```

## Generic HPC Usage

Write scheduler headers according to your cluster, activate the software
environment, and put the same `skmer-smk2 run` command in the script body. For
example, after your scheduler headers and environment activation:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 \
  --bowtie2-threads 2 \
  --fastp-threads 4 \
  --repair-threads 2 \
  --bbmerge-threads 2 \
  --total-mem-mb 180000 \
  --bowtie2-mem-mb 6000
```

Submit the script with your cluster's scheduler command. The workflow itself is
unchanged between local and HPC runs.

## Output Directory

All outputs are written under:

```text
results/
```

The most important final files are:

```text
results/stats/head_summary.sorted.tsv
results/stats/head_cutoff_candidates.tsv
results/skmer/tree.direct.tre
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
results/waster/waster.tree
results/waster/waster.branchlength.tree
results/mash/tree.direct.tre
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
results/mash/distance_heatmap.svg
```

When `-ref REF_FASTA` is used, the workflow also writes the final report/tree
files with a `.ref` suffix, for example:

```text
results/stats/head_summary.sorted.ref.tsv
results/skmer/tree.merged.ref.tre
results/waster/waster.branchlength.ref.tree
results/mash/tree.merged.ref.tre
results/mash/distance_heatmap.ref.svg
```

The unsuffixed files remain available as workflow intermediates and for
compatibility with earlier runs. The `.ref` files are the final outputs to use
when you want filenames to show that reference/plastid filtering was applied.

When `-skmer`, `-waster`, or `-mash` is used, only the selected analysis outputs
are required by the final workflow target. Shared preprocessing and statistics
outputs are still generated.

## Output Details

### Quality-Control Outputs

```text
results/<sample>/clean/<sample>.R1.clean.fq.gz
results/<sample>/clean/<sample>.R2.clean.fq.gz
results/<sample>/clean/<sample>.fastp.html
results/<sample>/clean/<sample>.fastp.json
```

These files come from `fastp`. The FASTQ files are cleaned reads. The HTML and
JSON reports summarize quality profiles, filtering rates, adapter trimming, and
read length distributions.

### Reference-Filtering Outputs

These files are produced only when `-ref` is supplied:

```text
results/ref/refDNA.*.bt2
results/<sample>/nDNA/un-conc-mate.1.fq.gz
results/<sample>/nDNA/un-conc-mate.2.fq.gz
results/<sample>/nDNA/<sample>.sam
results/<sample>/nDNA/un-conc-mate.1.fixed.fq.gz
results/<sample>/nDNA/un-conc-mate.2.fixed.fq.gz
results/<sample>/nDNA/singletons.fq.gz
```

`refDNA.*.bt2` is the Bowtie2 index. `un-conc-mate.*.fq.gz` contains read pairs
that did not align concordantly to the reference. The SAM file records
alignments to the reference. The `.fixed.fq.gz` files are synchronized pairs
after `repair.sh`; `singletons.fq.gz` contains reads whose mate was not retained.

### Merged And Concatenated Reads

```text
results/<sample>/sample_merged.fq
results/<sample>/sample_unmerged1.fq
results/<sample>/sample_unmerged2.fq
results/<sample>/sample.fq
results/<sample>/bbmerge.log
```

`sample_merged.fq` contains overlapping paired reads merged by `bbmerge.sh`.
`sample_unmerged1.fq` and `sample_unmerged2.fq` contain reads that could not be
merged. `sample.fq` concatenates merged and unmerged reads and is the input for
statistics and base-aware normalization.

If `bbmerge.sh` fails for a sample because of a BBMap/Java runtime error, does
not create all expected outputs, or exceeds `--bbmerge-timeout`, the workflow
falls back to a conservative unmerged-read mode for that sample:
`sample_merged.fq` is empty, and the cleaned R1/R2 reads are copied into
`sample_unmerged1.fq` and `sample_unmerged2.fq`. This keeps all reads available
for downstream k-mer analyses, but that sample will not benefit from overlap
merging. Check `bbmerge.log` for `BBMERGE_SUCCESS` or
`BBMERGE_FALLBACK_USED`.

### Statistics And Normalized FASTQ

```text
results/<sample>/stats/<sample>.post_filter.tsv
results/stats/post_filter_summary.sorted.tsv
results/stats/head_base_cutoff.txt
results/stats/head_cutoff_candidates.tsv
results/<sample>/stats/<sample>.head.tsv
results/stats/head_summary.sorted.tsv
results/<sample>/nDNAOK/<sample>.fq
```

`post_filter_summary.sorted.tsv` contains the per-sample read count, total bases,
and average read length before final head normalization. It is sorted from high
to low by total bases.

`head_base_cutoff.txt` records the base cutoff selected by `-s`.

`head_cutoff_candidates.tsv` compares multiple possible `-s` values and helps
decide whether a dataset is better normalized with 50, 75, 90, or another
percentile.

`head_summary.sorted.tsv` contains the same statistics after writing the final
normalized FASTQ files.

`nDNAOK/<sample>.fq` is the final normalized FASTQ used by Skmer, WASTER, and
Mash.

### Skmer Outputs

```text
results/skmer/dimtrx_main.txt
results/skmer/dimtrx_main_cor_.txt
results/skmer/dimtrx_main_cor_OK.phy
results/skmer/logs/reference.log
results/skmer/logs/subsample.log
results/skmer/logs/correct.log
results/skmer/bootstrap.trees
results/skmer/tree.direct.tre
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
```

`dimtrx_main.txt` is the initial Skmer distance matrix. `dimtrx_main_cor_.txt`
is the corrected Skmer distance matrix. `dimtrx_main_cor_OK.phy` is the PHYLIP
matrix passed to FastME.

`logs/reference.log`, `logs/subsample.log`, and `logs/correct.log` capture
Skmer's detailed output for the three Skmer calculation stages.

`tree.direct.tre` is the direct Skmer/FastME tree from the corrected distance
matrix.

`bootstrap.trees` contains all Skmer replicate trees. `tree.bootstrap.tre` is
the majority-rule consensus tree produced from those replicate trees. RAxML
support annotations are normalized from `:branch[support]` to standard internal
node labels so common tree viewers do not treat the support value as part of the
branch length.

`tree.merged.tre` is one Newick tree: the direct Skmer topology with matching
bootstrap support values copied onto internal nodes. This is the main Skmer tree
to inspect when you want direct branch lengths and bootstrap support together.

### WASTER Outputs

```text
results/waster/input.tsv
results/waster/waster.tree
results/waster/waster.branchlength.tree
results/waster/waster_branchlength.log
```

`input.tsv` is a two-column sample list generated by the workflow:

```text
sample_name    normalized_fastq_path
```

`waster.tree` is the WASTER topology tree inferred from the normalized FASTQ
files.

`waster.branchlength.tree` is produced by `waster_branchlength` using
`waster.tree` as the fixed topology and the same input FASTQ list to estimate
branch lengths. This is the WASTER tree to inspect when branch lengths are
needed. `waster_branchlength.log` records the branch-length estimation command
output.

### Mash Outputs

```text
results/mash/sketches/*.msh
results/mash/all.msh
results/mash/distances.tsv
results/mash/distances.phy
results/mash/distance_heatmap.svg
results/mash/bootstrap.trees
results/mash/tree.direct.tre
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
```

`sketches/*.msh` are per-sample Mash sketches. `all.msh` is the merged sketch
database. `distances.tsv` contains pairwise Mash distances. `distances.phy` is
the PHYLIP matrix passed to FastME.

`distance_heatmap.svg` visualizes the Mash distance matrix. It is useful for
checking outliers, unexpected sample similarity, or obvious sample mix-ups.

`tree.direct.tre` is the direct Mash/FastME tree. `bootstrap.trees` contains
replicate Mash trees generated from repeated sketches. `tree.bootstrap.tre` is
the normalized consensus tree. `tree.merged.tre` is one Newick tree with the
direct Mash topology and matching bootstrap support values on internal nodes.
This is the main Mash tree to inspect.

Older `skmer-smk2` runs may contain Mash tree files named
`mash_tree.direct.tre`, `mash_tree.bootstrap.tre`, and `mash_tree.merged.tre`.
Current versions automatically copy those legacy names to the standard
`tree.*.tre` names when reusing an existing work directory.

## Recommended Files To Inspect First

Quality and depth:

```text
input_qc/input_sample_report.tsv
input_qc/post_repair_sample_report.tsv
results/<sample>/clean/<sample>.fastp.html
results/stats/post_filter_summary.sorted.tsv
results/stats/head_cutoff_candidates.tsv
results/stats/head_summary.sorted.tsv
results/mash/distance_heatmap.svg
```

Final trees:

```text
results/skmer/tree.merged.tre
results/waster/waster.branchlength.tree
results/mash/tree.merged.tre
```

## Common Issues

### No Paired FASTQ Files Found

Make sure `-i` points directly to the FASTQ directory and that file names use
one of the supported R1/R2 naming styles.

### Reference Filtering Does Not Run

Bowtie2 filtering only runs when `-ref REF_FASTA` is provided. If `-ref` is
omitted, the workflow intentionally skips the reference-removal branch.

### PULP_CBC_CMD Or ILP Solver Warning

Older runs may print:

```text
Failed to solve scheduling problem with ILP solver, falling back to greedy scheduler
PULP_CBC_CMD: Not Available
```

This is a Snakemake scheduling warning, not a biological analysis error. It
means the optional CBC/ILP solver is unavailable, so Snakemake falls back to the
greedy scheduler. Current `skmer-smk2` uses `--scheduler greedy` by default to
avoid this warning. Use `--scheduler ilp` only if your environment has a working
CBC solver.

### Missing External Tool

Run:

```bash
skmer-smk2 doctor
```

Then install the missing tool or make sure it is available in `PATH`.

If the Skmer branch fails with `FileNotFoundError: ... 'seqtk'` or
`FileNotFoundError: ... 'jellyfish'`, install the missing program into the
active environment or make sure it is visible in `PATH`. Skmer calls these
programs internally during reference-distance estimation.

### Skmer Reference Uses Too Much Memory

If `rule skmer_reference` exits after a long run and the scheduler reports very
high memory use, lower the Skmer sketch size and optionally reduce Skmer's
internal thread count:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 \
  --skmer-sketch-size 50000 --skmer-threads 8 --skmer-mem-mb 180000
```

The default is `--skmer-sketch-size 100000 --skmer-threads 16`, matching the
original workflow. Smaller sketch sizes reduce memory and runtime pressure but
may slightly reduce Skmer distance resolution. Inspect
`results/skmer/logs/reference.log` for the detailed Skmer message.

### WASTER Is Killed With Exit Status 137

Exit status `137` usually means the scheduler or operating system killed WASTER
because memory was exhausted. For large datasets, run WASTER with fewer threads
and a larger memory scheduling estimate:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -waster \
  --waster-threads 2 \
  --waster-mem-mb 220000 \
  --total-mem-mb 240000
```

If the node still kills the job, request more memory from the scheduler, reduce
`--waster-threads` to `1`, or run only `-skmer -mash` first and run WASTER later
in a larger-memory job. WASTER is independent of the Skmer and Mash final trees.

### fastp Fails For One Sample

If the workflow stops at `rule fastp`, inspect the per-sample log:

```bash
cat results/<sample>/clean/<sample>.fastp.log
```

Common causes are truncated gzip files, incomplete FASTQ records, mismatched
R1/R2 files, or reads shorter than the configured filtering length. Check gzip
streams first:

```bash
gzip -t sample_1.fq.gz
gzip -t sample_2.fq.gz
```

For repairable FASTQ record problems, use the repair helper before rerunning the
workflow.

### bbmerge.sh Runtime Errors

If `bbmerge.sh` reports Java runtime errors such as `AssertionError` or
`IndexOutOfBoundsException`, or if one BBMerge job appears to run forever,
update to a recent `skmer-smk2` release and rerun the same command. The workflow
uses a wrapper that first tries BBMerge normally with `overwrite=t` and `-da`.
If BBMerge fails, logs a Java exception, misses expected outputs, or exceeds
`--bbmerge-timeout`, the wrapper writes a safe fallback: an empty merged-read
file plus the cleaned R1/R2 reads as unmerged outputs.

The default timeout is 14400 seconds per sample. Increase it for unusually large
datasets, for example:

```bash
skmer-smk2 run -i input_qc/input_for_skmer -ref /path/to/ref.fasta -s 75 -j 48 \
  --bbmerge-timeout 21600
```

Use `--bbmerge-timeout 0` only if you want to disable the timeout completely.

### Old Workflow Cache Is Used

The packaged Snakefile is copied into the working directory when the workflow
runs:

```text
.skmer_smk2_workflow/
```

After upgrading `skmer-smk2`, remove this directory if old commands still appear
in `--printshellcmds`.

### Existing Results Are Reused

Snakemake reuses completed outputs. For a clean rerun, use a new `--workdir` or
remove old `results/` after confirming that you do not need the previous output.
