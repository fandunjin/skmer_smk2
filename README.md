# skmer_smk2

`skmer_smk2` is a packaged Snakemake workflow for paired-end FASTQ processing, optional plastid-genome read removal, base-aware subsampling, and phylogenetic tree inference with Skmer, WASTER, and Mash.

The project is packaged as a Python command-line tool. Local runs and HPC runs use the same workflow and the same command; an HPC scheduler script is only an optional submission wrapper.

Demo reference and FASTQ files are copied from [`fandunjin/skmer_smk`](https://github.com/fandunjin/skmer_smk).

## What The Workflow Does

1. Detect paired FASTQ files named `_1/_2`, `_R1/_R2`, or `.R1/.R2`.
2. Run `fastp`.
3. Optionally remove plastid reads with `bowtie2` when `-ref` is provided.
4. Repair read pairs with `repair.sh`.
5. Merge reads with `bbmerge.sh`.
6. Summarize reads, bases, and average length before and after final subsampling.
7. Choose a base-count cutoff from the sorted sample-depth percentile.
8. Generate Skmer direct, bootstrap, and merged trees.
9. Generate a WASTER tree.
10. Generate Mash direct, bootstrap, merged trees, and a distance heatmap.

## Install

From the repository root:

```bash
python -m pip install .
```

For editable development:

```bash
python -m pip install -e .
```

The Python package installs the `skmer-smk2` command and bundles the Snakefile plus workflow helper scripts. It does not vendor large bioinformatics executables.

## Required External Tools

Make sure these tools are available in the active environment:

```text
snakemake
python
fastp
bowtie2
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

On HPC systems, install or load these with conda/modules, then run the same `skmer-smk2` command.

## Input Data

Put paired FASTQ files in one directory:

```text
SampleA_1.fq.gz
SampleA_2.fq.gz
SampleB_R1.fq.gz
SampleB_R2.fq.gz
SampleC.R1.fq.gz
SampleC.R2.fq.gz
```

For plastid removal, provide a FASTA reference:

```text
ref/refDNA.fasta
```

Bundled demo data:

```text
raw_data/ref.fna
raw_data/raw_data/sample*.R1.fq.gz
raw_data/raw_data/sample*.R2.fq.gz
```

## Run

Without plastid filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48
```

With plastid filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/refDNA.fasta -s 75 -j 48
```

Demo dry-run:

```bash
skmer-smk2 run -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 4 -b 2 --dry-run
```

The value of `-s` is the sorted sample-depth percentile used to choose the base-count cutoff. For example, `-s 75` sorts samples by total bases from large to small, takes the base count at the 75% position, and uses that value for final FASTQ truncation.

## Optional Templates

Write optional helper templates:

```bash
skmer-smk2 init -o run_templates
```

This writes:

```text
run_templates/jsub_submit.sh
run_templates/scan_repair_fastq.sh
```

To also export the packaged Snakefile and helper scripts:

```bash
skmer-smk2 init -o run_templates --with-workflow
```

## HPC Usage

HPC uses the same packaged workflow. The only difference is that the scheduler submits a wrapper script.

Example using the bundled JSUB template:

```bash
skmer-smk2 init -o run_templates

export WORKDIR=/path/to/workdir
export INPUT_DIR=/path/to/fastq_dir
export REF=/path/to/refDNA.fasta
export CONDA_PROFILE=/path/to/conda/etc/profile.d/conda.sh
export CONDA_ENV=your_bioinfo_env
export THREADS=48
export SAMPLE_PERCENTILE=75
export BOOTSTRAPS=100
export EXCLUDE_SAMPLES=""

cd "${WORKDIR}"
jsub < /path/to/run_templates/jsub_submit.sh
```

For schedulers other than JSUB, adapt only the scheduler header and submission command. The inner workflow command should remain `skmer-smk2 run`.

## FASTQ Scan And Repair

Copy the repair helper:

```bash
skmer-smk2 repair-fastq --workdir . --copy-only
```

Or run it directly if `bash` and `seqkit` are available:

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
