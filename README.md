# skmer_smk2

`skmer_smk2` packages a Snakemake workflow for paired-end FASTQ processing, optional plastid-genome read removal, base-aware subsampling, and tree inference with Skmer, WASTER, and Mash.

Local runs and HPC runs use the same command:

```bash
skmer-smk2 run -i FASTQ_DIR -ref REF_FASTA -s 75 -j 48
```

HPC is not a separate workflow. The scheduler only launches this same command on a compute node.

Demo reference and FASTQ files are copied from [`fandunjin/skmer_smk`](https://github.com/fandunjin/skmer_smk).

## Install

Install from GitHub:

```bash
python -m pip install git+https://github.com/fandunjin/skmer_smk2.git
```

Or install from a local clone:

```bash
git clone https://github.com/fandunjin/skmer_smk2.git
cd skmer_smk2
python -m pip install .
```

For editable development:

```bash
python -m pip install -e .
```

The Python package supports Python 3.8 or newer. It installs the `skmer-smk2` command and bundles the Snakefile plus workflow helper scripts. It does not force-install Snakemake or other large bioinformatics executables, so it can be installed into existing HPC environments without changing their software stack.

## Check Environment

Check the active environment:

```bash
skmer-smk2 doctor
```

`doctor` is advisory by default: it prints OK/WARN rows but does not fail just because a tool was not found by the current shell. This is useful on HPC systems where modules or job scripts may initialize PATH differently. Use strict mode only when you want missing tools to return a non-zero exit code:

```bash
skmer-smk2 doctor --strict
```

Show help:

```bash
skmer-smk2 -h
skmer-smk2 run -h
skmer-smk2 doctor -h
```

`doctor` checks:

```text
snakemake python fastp bowtie2 repair.sh bbmerge.sh skmer fastme raxmlHPC waster mash seqkit gzip
```

Install missing conda-available packages into the current environment:

```bash
skmer-smk2 doctor --install
```

`doctor --install` prefers `mamba` and falls back to `conda`. It uses `conda-forge` and `bioconda`, prints the packages before installing, and marks tools that need manual installation. Existing environments that already contain the required tools can skip this step.

## Run

FASTQ files can be named with `_1/_2`, `_R1/_R2`, or `.R1/.R2`:

```text
SampleA_1.fq.gz
SampleA_2.fq.gz
SampleB_R1.fq.gz
SampleB_R2.fq.gz
SampleC.R1.fq.gz
SampleC.R2.fq.gz
```

Run without plastid filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48
```

Run with plastid filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/refDNA.fasta -s 75 -j 48
```

Demo dry-run:

```bash
skmer-smk2 run -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 1 -b 2 --dry-run
```

The value of `-s` is the sorted sample-depth percentile used to choose the base-count cutoff. For example, `-s 75` sorts samples by total bases from large to small, takes the base count at the 75% position, and uses that value for final FASTQ truncation.

## HPC Example

Write the scheduler headers required by your cluster, activate the environment, then run the same command:

```bash
#!/bin/bash
# scheduler headers written by user

source /path/to/conda.sh
conda activate your_env

skmer-smk2 run -i /path/to/fq -ref /path/to/ref.fa -s 75 -j 48
```

Submit the script with your cluster command, for example `jsub < skmer.sh`, `sbatch skmer.sh`, or `qsub skmer.sh`.

## FASTQ Scan And Repair

Copy the repair helper:

```bash
skmer-smk2 repair-fastq --workdir . --copy-only
```

Run it directly if `bash` and `seqkit` are available:

```bash
skmer-smk2 repair-fastq -i /path/to/fastq_dir --workdir .
```

The report is written under:

```text
repaired_fastq/fastq_repair_report.tsv
repaired_fastq/input_for_skmer/
```

Report actions:

```text
UNCHANGED       valid input; no repaired copy kept
REPAIRED        repaired by seqkit sana
NEED_REUPLOAD   gzip stream is truncated or physically incomplete
FAILED          seqkit sana failed
```

Files marked `NEED_REUPLOAD` should be uploaded again from the original source.

## Main Outputs

```text
results/stats/post_filter_summary.sorted.tsv
results/stats/head_summary.sorted.tsv
results/skmer/tree.direct.tre
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
results/waster/waster.tree
results/mash/tree.direct.tre
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
results/mash/distance_heatmap.svg
```

## Optional Helpers

Write optional helper scripts:

```bash
skmer-smk2 init -o run_templates
```

This writes:

```text
run_templates/jsub_submit.sh
run_templates/scan_repair_fastq.sh
```

These helpers are examples only. The main interface remains `skmer-smk2 run ...`.

## Repository Layout

```text
.
|-- pyproject.toml
|-- skmer_smk2/
|   |-- cli.py
|   |-- templates/
|   |   |-- jsub_submit.sh
|   |   `-- scan_repair_fastq.sh
|   `-- workflow/
|       |-- Snakefile
|       `-- scripts/
|-- raw_data/
`-- README.md
```
