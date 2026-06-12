# skmer-smk2

`skmer-smk2` is a packaged Snakemake workflow for phylogenetic analysis from paired-end FASTQ reads. It combines read quality control, optional reference-based read removal, read merging, base-aware downsampling, and tree inference with Skmer, WASTER, and Mash behind one command-line interface:

```bash
skmer-smk2 run -i FASTQ_DIR -ref REF_FASTA -s 75 -j 48
```

The package is designed for both local Linux workstations and HPC clusters. The workflow command is the same in both cases; an HPC scheduler script only needs to request resources, activate the software environment, and launch `skmer-smk2 run`.

## What The Workflow Does

Starting from paired-end FASTQ files, `skmer-smk2` performs the following operations:

1. Discovers complete R1/R2 FASTQ pairs in the input directory.
2. Runs `fastp` for adapter trimming, quality filtering, and paired-end correction.
3. Optionally removes reads that map to a reference genome with `bowtie2`. This is useful for excluding plastid, chloroplast, mitochondrial, or other reference-derived reads.
4. Repairs pair relationships after reference filtering with BBMap `repair.sh`.
5. Merges overlapping paired reads with BBMap `bbmerge.sh`.
6. Concatenates merged and unmerged reads into one per-sample FASTQ.
7. Computes per-sample read and base statistics.
8. Selects a base-count cutoff from a user-defined sample-depth percentile.
9. Truncates each sample to the same base-count scale.
10. Runs Skmer, WASTER, and Mash phylogenetic analyses.
11. Builds direct trees, bootstrap trees, merged consensus trees, and a Mash distance heatmap.

## Repository Layout

```text
.
|-- README.md
|-- pyproject.toml
|-- src/
|   `-- skmer_smk2/
|       |-- cli.py
|       |-- templates/
|       |   |-- jsub_submit.sh
|       |   `-- scan_repair_fastq.sh
|       `-- workflow/
|           |-- Snakefile
|           `-- scripts/
|-- raw_data/
|-- repair_head_truncated_fastq.sh
`-- replace_repaired_fastq.sh
```

Important files and directories:

```text
src/skmer_smk2/cli.py
```

The command-line interface. It provides `run`, `doctor`, `init`, and `repair-fastq`.

```text
src/skmer_smk2/workflow/Snakefile
```

The packaged Snakemake workflow used by `skmer-smk2 run`.

```text
src/skmer_smk2/workflow/scripts/
```

Helper scripts for statistics, downsampling, distance matrix conversion, tree merging, and Mash heatmap plotting.

```text
raw_data/
```

Small bundled demo files for testing installation and dry runs.

```text
repair_head_truncated_fastq.sh
replace_repaired_fastq.sh
```

Standalone helper scripts for repairing problematic FASTQ files and safely replacing confirmed repaired files.

## Requirements

The Python package installs the `skmer-smk2` command and the bundled workflow files. Large bioinformatics tools are intentionally kept as external dependencies so the package can be installed into existing HPC environments without modifying the cluster software stack.

Required or commonly used external tools:

```text
python >= 3.8
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
mash
seqkit
gzip
```

Most dependencies can be installed from `conda-forge` and `bioconda`. `skmer` and `waster` may need manual installation from their upstream projects, depending on your environment.

## Installation

### Install From GitHub

```bash
python -m pip install git+https://github.com/fandunjin/skmer_smk2.git
```

If `git clone` over HTTPS is unstable on an HPC login node, install from the GitHub zip archive instead:

```bash
python -m pip install --no-cache-dir --force-reinstall \
  https://github.com/fandunjin/skmer_smk2/archive/refs/heads/main.zip
```

### Install From A Local Clone

```bash
git clone https://github.com/fandunjin/skmer_smk2.git
cd skmer_smk2
python -m pip install .
```

For development:

```bash
python -m pip install -e .
```

Confirm the command is available:

```bash
skmer-smk2 --version
skmer-smk2 -h
```

### Install External Tools With Conda Or Mamba

Activate the environment you want to use:

```bash
conda activate your_env
```

Install conda-available tools:

```bash
mamba install -c conda-forge -c bioconda \
  snakemake fastp bowtie2 bbmap fastme raxml mash seqkit gzip
```

If `mamba` is unavailable, use `conda`:

```bash
conda install -c conda-forge -c bioconda \
  snakemake fastp bowtie2 bbmap fastme raxml mash seqkit gzip
```

Then install `skmer` and `waster` following their upstream instructions, and make sure both executables are in `PATH`:

```bash
which skmer
which waster
```

## Environment Check

Run:

```bash
skmer-smk2 doctor
```

`doctor` checks whether the required tools can be resolved from the current shell. By default it is advisory: it prints `OK` or `WARN` rows but does not fail only because some tools are missing. This is useful on clusters where scheduler scripts may initialize `PATH` differently from login shells.

Use strict mode when you want missing tools to return a non-zero exit code:

```bash
skmer-smk2 doctor --strict
```

Install missing conda-available packages into the active environment:

```bash
skmer-smk2 doctor --install
```

Tools checked by `doctor`:

```text
snakemake python fastp bowtie2 bowtie2-build repair.sh bbmerge.sh
skmer fastme raxmlHPC waster mash seqkit gzip
```

## Input FASTQ Naming

The input directory must directly contain paired FASTQ files. Do not pass the repository root or a parent directory.

Supported paired-end naming styles:

```text
SampleA_1.fq.gz        SampleA_2.fq.gz
SampleB_R1.fq.gz       SampleB_R2.fq.gz
SampleC.R1.fastq.gz    SampleC.R2.fastq.gz
SampleD-R1.fq.gz       SampleD-R2.fq.gz
```

The sample name is the shared prefix before the mate suffix. For example:

```text
H_asiatica_SAMC1020836_1.fq.gz
H_asiatica_SAMC1020836_2.fq.gz
```

is discovered as sample:

```text
H_asiatica_SAMC1020836
```

Before running the workflow, confirm that every sample has both mates:

```bash
ls *_1.fq.gz *_2.fq.gz
```

## Basic Usage

Run without reference filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 --printshellcmds
```

Run with reference filtering:

```bash
skmer-smk2 run \
  -i /path/to/fastq_dir \
  -ref /path/to/refDNA.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

Dry run:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/refDNA.fasta -s 75 -j 1 --dry-run
```

Pass extra arguments to Snakemake after `--`:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 -- --keep-going
```

### Main Run Options

```text
-i, --input
```

Directory containing paired FASTQ files.

```text
-ref, --ref
```

Optional reference FASTA for `bowtie2` filtering. Reads mapping to this reference are removed from downstream analysis.

```text
-s, --sample-percentile
```

Sorted sample-depth percentile used to choose the base-count cutoff. Samples are sorted by total post-filter bases from high to low; the selected percentile determines how many bases each sample contributes to final analysis. The default is `75`.

```text
-j, --jobs
```

Number of cores/jobs available to Snakemake.

```text
-b, --bootstraps
```

Number of bootstrap replicates. Default: `100`.

```text
--exclude-samples
```

Comma- or whitespace-separated sample names to skip.

```text
--workdir
```

Run directory for `results/` and the materialized workflow cache. Default: current directory.

```text
--latency-wait
```

Snakemake latency wait in seconds. Default: `120`.

```text
--dry-run
```

Preview the workflow without running commands.

```text
--printshellcmds
```

Print all shell commands executed by Snakemake.

## Demo Data

The repository contains small demo FASTQ and reference files:

```text
raw_data/raw_data/
raw_data/ref.fna
```

Run a demo dry run:

```bash
skmer-smk2 run -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 1 -b 2 --dry-run
```

## HPC Usage

On an HPC cluster, write a scheduler script that activates the environment and then runs the same workflow command.

Example `jsub` script:

```bash
#!/bin/bash
#JSUB -J skmer
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_log.%J
#JSUB -e skmer_err.%J
#JSUB -cwd .

source /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh
conda activate 01bio

skmer-smk2 run \
  -i /hpcfile/users/92024286/Huperzia \
  -ref /path/to/refDNA.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

Submit:

```bash
jsub < skmer_hpc.sh
```

For `sbatch`, `qsub`, or another scheduler, replace only the scheduler header. The `skmer-smk2 run` command remains the same.

## FASTQ Scan And Repair

FASTQ corruption, truncated records, and mismatched R1/R2 files should be detected before the main workflow. `skmer-smk2` provides a general repair helper, and this repository also includes a specialized helper for head/truncated FASTQ cases.

### General FASTQ Scan

Copy the repair helper:

```bash
skmer-smk2 repair-fastq --workdir . --copy-only
```

Run it:

```bash
bash scan_repair_fastq.sh /path/to/fastq_dir
```

Outputs:

```text
repaired_fastq/fastq_repair_report.tsv
repaired_fastq/input_for_skmer/
```

Report actions:

```text
UNCHANGED
```

The original file passed checks. A symlink to the original file is placed in `input_for_skmer/`.

```text
REPAIRED
```

`seqkit sana` produced a repaired copy and the read count changed. The repaired file is linked into `input_for_skmer/`.

```text
NEED_REUPLOAD
```

The gzip stream is incomplete or corrupted. The safest fix is to download or upload the file again from the original source.

```text
FAILED
```

Repair failed. Check the generated log file.

### Pair-Aware Repair For Truncated FASTQ

The standalone script:

```bash
repair_head_truncated_fastq.sh
```

is intended for cases where individual FASTQ records are malformed but the gzip stream is still readable. It runs `seqkit sana` on both mates, then runs `seqkit pair`, and finally validates that the output is usable as paired-end data.

Run directly:

```bash
bash repair_head_truncated_fastq.sh
```

Or submit to a `jsub` cluster:

```bash
jsub < repair_head_truncated_fastq.sh
```

Outputs:

```text
fastq_repair_report.tsv
repaired_for_skmer/
failed_repair_outputs/
repair_tmp/
```

Status values:

```text
REPAIRED
```

The pair was repaired successfully. R1 and R2 are valid gzip files, both have more than zero reads, and read counts are equal.

```text
EMPTY_OUTPUT
```

`seqkit pair` found zero paired reads. This usually means the R1 and R2 read identifiers do not match, or the mate files are from different datasets.

```text
UNPAIRED_OUTPUT
```

R1 and R2 output counts differ. Do not use this pair as paired-end input without further investigation.

```text
BAD_OUTPUT
```

The repaired files failed gzip or FASTQ statistics checks.

```text
MISSING_INPUT
```

One or both source files were missing or empty.

### Replacing Original FASTQ With Confirmed Repairs

Use:

```bash
replace_repaired_fastq.sh
```

only for samples that were confirmed as `REPAIRED`. The script:

1. Checks that the repaired files exist.
2. Runs `gzip -t`.
3. Confirms R1 and R2 read counts are equal and greater than zero.
4. Creates a timestamped backup directory such as `original_fastq_backup_20260612_133000`.
5. Moves the original FASTQ files into the backup directory.
6. Copies repaired files from `repaired_for_skmer/` back to the original filenames.

Run:

```bash
bash replace_repaired_fastq.sh
```

Verify:

```bash
ls -lh original_fastq_backup_*
ls -lh *.fq.gz
```

## Output Directory Overview

The workflow writes results under:

```text
results/
```

High-level layout:

```text
results/
|-- <sample>/
|   |-- clean/
|   |-- nDNA/
|   |-- nDNAOK/
|   `-- stats/
|-- stats/
|-- skmer/
|-- skmer_input_dir/
|-- waster/
`-- mash/
```

## Detailed Output Interpretation

### Per-Sample Clean Reads

```text
results/<sample>/clean/<sample>.R1.clean.fq.gz
results/<sample>/clean/<sample>.R2.clean.fq.gz
results/<sample>/clean/<sample>.fastp.json
results/<sample>/clean/<sample>.fastp.html
```

These are the `fastp` outputs.

The cleaned FASTQ files are adapter-trimmed, quality-filtered paired reads. They are the input for optional reference filtering or, if no reference is supplied, for read merging.

The HTML and JSON reports summarize quality before and after filtering, adapter trimming, read length distributions, duplication, and filtering reasons. Use these reports to check whether a sample has unusually poor quality, too many short reads, or excessive filtering.

### Reference-Filtered Reads

These files are produced only when `-ref` is supplied:

```text
results/<sample>/nDNA/un-conc-mate.1.fq.gz
results/<sample>/nDNA/un-conc-mate.2.fq.gz
results/<sample>/nDNA/<sample>.sam
```

`bowtie2` writes reads that did not align concordantly to the reference as `un-conc-mate.*.fq.gz`. These are the reads retained for downstream phylogenetic analysis.

The SAM file records alignments to the reference. It can be inspected if you need to understand how many reads mapped to the reference or why a sample lost many reads.

### Pair-Repaired Reference-Filtered Reads

```text
results/<sample>/nDNA/un-conc-mate.1.fixed.fq.gz
results/<sample>/nDNA/un-conc-mate.2.fixed.fq.gz
results/<sample>/nDNA/singletons.fq.gz
```

After reference filtering, some mates may be removed independently. `repair.sh` restores synchronized paired FASTQ files. The `.fixed.fq.gz` files are used by `bbmerge.sh`.

`singletons.fq.gz` contains reads whose mate was missing after filtering. These reads are not used as paired reads in the downstream merge step, but the file is useful for diagnosing heavy mate loss.

### Merged And Unmerged Reads

```text
results/<sample>/sample_merged.fq
results/<sample>/sample_unmerged1.fq
results/<sample>/sample_unmerged2.fq
results/<sample>/sample.fq
```

`bbmerge.sh` attempts to merge overlapping R1/R2 reads.

`sample_merged.fq` contains merged read pairs. `sample_unmerged1.fq` and `sample_unmerged2.fq` contain reads that could not be merged. `sample.fq` concatenates all three files and becomes the per-sample read set for statistics and base-aware downsampling.

### Post-Filter Statistics

```text
results/<sample>/stats/<sample>.post_filter.tsv
results/stats/post_filter_summary.sorted.tsv
```

The per-sample `post_filter.tsv` files report read and base counts after cleaning, optional reference filtering, pair repair, merging, and concatenation.

`post_filter_summary.sorted.tsv` combines all samples and sorts them by total bases. This file is important for checking whether all samples have enough sequence data and whether some samples are much shallower than the others.

### Base Cutoff

```text
results/stats/head_base_cutoff.txt
```

This file stores the base-count cutoff selected from `post_filter_summary.sorted.tsv` using the `-s/--sample-percentile` value.

For example, `-s 75` sorts samples from high to low total bases and selects the base count at the 75% position. That cutoff is then used to truncate every sample to a comparable amount of sequence.

This step reduces bias caused by unequal sequencing depth.

### Final Downsampled FASTQ

```text
results/<sample>/nDNAOK/<sample>.fq
results/<sample>/stats/<sample>.head.tsv
results/stats/head_summary.sorted.tsv
```

`nDNAOK/<sample>.fq` is the final per-sample FASTQ used by Skmer, WASTER, and Mash.

`head.tsv` records the final read and base count for each sample after downsampling.

`head_summary.sorted.tsv` summarizes final input depth across all samples. This is one of the most important QC files. All retained samples should have comparable final base counts unless the sample had fewer bases than the selected cutoff.

### Skmer Input Directory

```text
results/skmer_input_dir/
```

This directory contains symlinks, hard links, or copies of the final `nDNAOK/<sample>.fq` files. It is created because Skmer expects an input directory containing one file per sample.

If a platform does not support symlinks, the workflow falls back to hard links or file copies.

### Skmer Distance And Trees

```text
results/skmer/dimtrx_main.txt
results/skmer/dimtrx_main_cor_.txt
results/skmer/dimtrx_main_cor_OK.phy
results/skmer/tree.direct.tre
results/skmer/bootstrap.trees
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
results/skmer/subsample/
```

Meaning:

```text
dimtrx_main.txt
```

The main Skmer distance matrix before correction.

```text
dimtrx_main_cor_.txt
```

The corrected Skmer distance matrix generated by `skmer correct`.

```text
dimtrx_main_cor_OK.phy
```

The corrected distance matrix converted to PHYLIP format for tree inference with `fastme`.

```text
tree.direct.tre
```

The direct FastME tree inferred from the corrected Skmer distance matrix.

```text
bootstrap.trees
```

All bootstrap replicate trees concatenated into one file.

```text
tree.bootstrap.tre
```

The RAxML majority-rule extended consensus tree from bootstrap replicate trees.

```text
tree.merged.tre
```

The direct Skmer tree annotated or merged with bootstrap support information. This is often the most useful Skmer tree for reporting.

```text
subsample/
```

Skmer bootstrap replicate distance matrices and trees.

### WASTER Output

```text
results/waster/input.tsv
results/waster/waster.tree
```

`input.tsv` lists sample names and final FASTQ paths passed to WASTER.

`waster.tree` is the WASTER-inferred tree. It provides an independent tree estimate from the same downsampled input reads. Compare it with Skmer and Mash trees to assess whether the major topology is robust across methods.

### Mash Distances, Heatmap, And Trees

```text
results/mash/sketches/
results/mash/all.msh
results/mash/distances.tsv
results/mash/distances.phy
results/mash/distance_heatmap.svg
results/mash/tree.direct.tre
results/mash/bootstrap.trees
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
results/mash/bootstrap/
```

Meaning:

```text
sketches/
```

Per-sample Mash sketches built from final downsampled FASTQ files.

```text
all.msh
```

Combined Mash sketch database.

```text
distances.tsv
```

Pairwise Mash distances among all samples. Lower values indicate more similar samples.

```text
distances.phy
```

Mash distance matrix converted to PHYLIP format for FastME.

```text
distance_heatmap.svg
```

Visual heatmap of the Mash distance matrix. This is useful for quickly detecting clustering, outliers, mislabeled samples, or unexpectedly distant samples.

```text
tree.direct.tre
```

FastME tree inferred directly from the Mash distance matrix.

```text
bootstrap.trees
```

Concatenated Mash bootstrap replicate trees.

```text
tree.bootstrap.tre
```

RAxML majority-rule extended consensus tree from Mash bootstrap trees.

```text
tree.merged.tre
```

Mash direct tree merged with bootstrap support information. This is often the most useful Mash tree for reporting.

```text
bootstrap/
```

Mash bootstrap replicate sketches, distance matrices, and trees.

## Which Outputs Should I Use?

For quality control:

```text
results/<sample>/clean/<sample>.fastp.html
results/stats/post_filter_summary.sorted.tsv
results/stats/head_summary.sorted.tsv
results/mash/distance_heatmap.svg
```

For final Skmer phylogeny:

```text
results/skmer/tree.merged.tre
```

For a Skmer tree without bootstrap annotation:

```text
results/skmer/tree.direct.tre
```

For the Skmer bootstrap consensus:

```text
results/skmer/tree.bootstrap.tre
```

For comparison with another method:

```text
results/waster/waster.tree
results/mash/tree.merged.tre
```

For distance-matrix based inspection:

```text
results/skmer/dimtrx_main_cor_OK.phy
results/mash/distances.tsv
results/mash/distance_heatmap.svg
```

## Common Problems

### `ADDR2LINE: unbound variable`

Some conda `activate.d` scripts reference unset variables. If a shell script enables `set -u` before `conda activate`, activation may fail with an error like:

```text
ADDR2LINE: unbound variable
```

Temporarily disable nounset during activation:

```bash
set +u
source /path/to/conda.sh
conda activate your_env
set -u
```

The helper scripts in this repository already apply this pattern.

### `No paired FASTQ files found`

Check that `-i` points to the directory containing FASTQ files and that filenames match a supported naming style:

```text
sample_1.fq.gz / sample_2.fq.gz
sample_R1.fq.gz / sample_R2.fq.gz
sample.R1.fq.gz / sample.R2.fq.gz
sample-R1.fq.gz / sample-R2.fq.gz
```

### `seqkit pair` Saves Zero Paired Reads

If `seqkit sana` reports many passing records but `seqkit pair` reports:

```text
0 paired-end reads saved
```

then R1 and R2 read identifiers probably do not overlap. This usually means:

- The R1 and R2 files are from different samples.
- One mate file was downloaded or uploaded incorrectly.
- The data source is not a true paired-end dataset.
- Read names were rewritten differently between mates.

This cannot be fixed by ordinary repair. Re-download or re-upload the correct mate files.

### Old Workflow Cache Is Still Used

`skmer-smk2 run` materializes the packaged workflow into:

```text
.skmer_smk2_workflow/
```

If you upgrade the package but still see old command paths in `--printshellcmds`, remove the cache:

```bash
rm -rf .skmer_smk2_workflow
```

Then rerun the workflow.

### Existing Results Are Reused

Snakemake reuses existing complete outputs. If you want a clean rerun, use a new `--workdir` or carefully remove previous `results/` and `.skmer_smk2_workflow/` directories after confirming the paths.

## Recommended End-To-End Workflow

1. Activate the environment:

```bash
conda activate your_env
```

2. Check dependencies:

```bash
skmer-smk2 doctor
```

3. Validate or repair FASTQ files:

```bash
bash repair_head_truncated_fastq.sh
cat fastq_repair_report.tsv
```

4. Keep only samples that are valid paired-end data.

5. Run the main workflow:

```bash
skmer-smk2 run \
  -i /path/to/fastq_dir \
  -ref /path/to/refDNA.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

6. Inspect QC outputs:

```bash
cat results/stats/post_filter_summary.sorted.tsv
cat results/stats/head_summary.sorted.tsv
```

7. Use final trees:

```text
results/skmer/tree.merged.tre
results/mash/tree.merged.tre
results/waster/waster.tree
```

