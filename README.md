# skmer-smk2

`skmer-smk2` is a packaged Snakemake workflow for building phylogenetic trees
from paired-end FASTQ reads. It provides one command for read cleaning, optional
reference-based plastid filtering, base-aware depth normalization, and tree
inference with Skmer, WASTER, and Mash.

```bash
skmer-smk2 run -i FASTQ_DIR -ref REF_FASTA -s 75 -j 48
```

The same command can be used on a local workstation or inside an HPC scheduler
script. HPC does not use a different workflow; the scheduler only submits this
command to a compute node.

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
| WASTER tree | one WASTER result tree |
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
mash
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
`skmer` and `waster` may need manual installation depending on your environment.

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

## Running The Workflow

Run without reference filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48
```

Run with reference filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48
```

Run only selected analysis branches:

```bash
# Mash only
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -mash

# Skmer and Mash, without WASTER
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -skmer -mash

# WASTER only
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48 -waster
```

If none of `-skmer`, `-waster`, or `-mash` is supplied, all three branches run.

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
| `-skmer` | Run the Skmer branch |
| `-waster` | Run the WASTER branch |
| `-mash` | Run the Mash branch |
| `-j`, `--jobs` | Snakemake cores/jobs |
| `-b`, `--bootstraps` | Bootstrap replicate count for Skmer and Mash; default `100` |
| `--exclude-samples` | Comma- or whitespace-separated sample names to skip |
| `--workdir` | Directory where `results/` and workflow cache are written; default current directory |
| `--latency-wait` | Snakemake latency wait seconds; default `120` |
| `--dry-run` | Build and print the DAG without running jobs |
| `--printshellcmds` | Print shell commands from Snakemake |

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
environment, and put the same `skmer-smk2 run` command in the script body.

```bash
#!/bin/bash
# scheduler headers go here

source /path/to/conda.sh
conda activate your_environment

skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/ref.fasta -s 75 -j 48
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
results/mash/tree.direct.tre
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
results/mash/distance_heatmap.svg
```

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
```

`sample_merged.fq` contains overlapping paired reads merged by `bbmerge.sh`.
`sample_unmerged1.fq` and `sample_unmerged2.fq` contain reads that could not be
merged. `sample.fq` concatenates merged and unmerged reads and is the input for
statistics and base-aware normalization.

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
results/skmer/bootstrap.trees
results/skmer/tree.direct.tre
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
```

`dimtrx_main.txt` is the initial Skmer distance matrix. `dimtrx_main_cor_.txt`
is the corrected Skmer distance matrix. `dimtrx_main_cor_OK.phy` is the PHYLIP
matrix passed to FastME.

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
```

`input.tsv` is a two-column sample list generated by the workflow:

```text
sample_name    normalized_fastq_path
```

`waster.tree` is the WASTER tree inferred from the normalized FASTQ files.

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

## Recommended Files To Inspect First

Quality and depth:

```text
results/<sample>/clean/<sample>.fastp.html
results/stats/post_filter_summary.sorted.tsv
results/stats/head_cutoff_candidates.tsv
results/stats/head_summary.sorted.tsv
results/mash/distance_heatmap.svg
```

Final trees:

```text
results/skmer/tree.merged.tre
results/waster/waster.tree
results/mash/tree.merged.tre
```

## FASTQ Repair Helper

The package also includes a helper for scanning and repairing FASTQ files before
running the main workflow.

Copy the helper script:

```bash
skmer-smk2 repair-fastq --workdir . --copy-only
```

Run it on a FASTQ directory:

```bash
bash scan_repair_fastq.sh /path/to/fastq_dir
```

Main outputs:

```text
repaired_fastq/fastq_repair_report.tsv
repaired_fastq/input_for_skmer/
```

Report actions:

| Action | Meaning |
| --- | --- |
| `UNCHANGED` | The original file passed checks |
| `REPAIRED` | A repaired copy was produced by `seqkit sana` |
| `NEED_REUPLOAD` | The gzip stream is incomplete or physically truncated |
| `FAILED` | Repair failed; inspect the log |

Use `repaired_fastq/input_for_skmer/` as the `-i` directory when it contains the
files you want to analyze.

## Common Issues

### No Paired FASTQ Files Found

Make sure `-i` points directly to the FASTQ directory and that file names use
one of the supported R1/R2 naming styles.

### Reference Filtering Does Not Run

Bowtie2 filtering only runs when `-ref REF_FASTA` is provided. If `-ref` is
omitted, the workflow intentionally skips the reference-removal branch.

### Missing External Tool

Run:

```bash
skmer-smk2 doctor
```

Then install the missing tool or make sure it is available in `PATH`.

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
