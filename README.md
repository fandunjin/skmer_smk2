# skmer-smk2

`skmer-smk2` is a packaged Snakemake workflow for phylogenetic analysis from paired-end FASTQ reads. It wraps read cleaning, optional reference-based filtering, depth normalization, and tree inference with Skmer, WASTER, and Mash.

```bash
skmer-smk2 run -i FASTQ_DIR -ref REF_FASTA -s 75 -j 48
```

## Workflow

For each paired-end sample, the workflow:

1. Detects complete R1/R2 FASTQ pairs.
2. Cleans reads with `fastp`.
3. Optionally removes reads mapping to a reference genome with `bowtie2`.
4. Repairs paired reads with BBMap `repair.sh`.
5. Merges overlapping reads with BBMap `bbmerge.sh`.
6. Computes per-sample read and base statistics.
7. Normalizes samples to a shared base-count cutoff.
8. Builds phylogenetic outputs with Skmer, WASTER, and Mash.

## Installation

Install the Python package:

```bash
python -m pip install git+https://github.com/fandunjin/skmer_smk2.git
```

Or install from a local clone:

```bash
git clone https://github.com/fandunjin/skmer_smk2.git
cd skmer_smk2
python -m pip install .
```

For development:

```bash
python -m pip install -e .
```

Check the installed command:

```bash
skmer-smk2 --version
skmer-smk2 -h
```

## Dependencies

The package installs the command-line wrapper and bundled Snakemake workflow. External bioinformatics tools must be available in the active environment:

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

Most tools can be installed with conda or mamba:

```bash
mamba install -c conda-forge -c bioconda \
  snakemake fastp bowtie2 bbmap fastme raxml mash seqkit gzip
```

Install `skmer` and `waster` according to their upstream instructions if they are not available in your package manager.

Check the active environment:

```bash
skmer-smk2 doctor
```

Use strict mode when missing tools should return a non-zero exit code:

```bash
skmer-smk2 doctor --strict
```

Install missing conda-available tools into the active environment:

```bash
skmer-smk2 doctor --install
```

## Input

The input directory must directly contain paired FASTQ files. Supported naming styles:

```text
SampleA_1.fq.gz      SampleA_2.fq.gz
SampleB_R1.fq.gz     SampleB_R2.fq.gz
SampleC.R1.fq.gz     SampleC.R2.fq.gz
SampleD-R1.fq.gz     SampleD-R2.fq.gz
```

The sample name is the shared prefix before the mate suffix. For example, `SampleA_1.fq.gz` and `SampleA_2.fq.gz` are detected as sample `SampleA`.

## Usage

Run without reference filtering:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48
```

Run with reference filtering:

```bash
skmer-smk2 run \
  -i /path/to/fastq_dir \
  -ref /path/to/reference.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

Dry run:

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/reference.fasta -s 75 -j 1 --dry-run
```

Useful options:

```text
-i, --input              Directory containing paired FASTQ files.
-ref, --ref              Optional reference FASTA for read removal.
-s, --sample-percentile  Percentile used to choose the base-count cutoff. Default: 75.
-j, --jobs               Number of Snakemake cores/jobs.
-b, --bootstraps         Number of bootstrap replicates. Default: 100.
--exclude-samples        Comma- or whitespace-separated sample names to skip.
--workdir                Working directory for results. Default: current directory.
--dry-run                Preview the workflow.
--printshellcmds         Print commands executed by Snakemake.
```

Extra Snakemake arguments can be passed after `--`:

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 -- --keep-going
```

## Example Data

The repository includes small test data:

```text
raw_data/raw_data/
raw_data/ref.fna
```

Run a dry run:

```bash
skmer-smk2 run -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 1 -b 2 --dry-run
```

## FASTQ Repair Helper

Copy the repair helper script:

```bash
skmer-smk2 repair-fastq --workdir . --copy-only
```

Run it on a FASTQ directory:

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
UNCHANGED      The original file passed checks.
REPAIRED       A repaired copy was produced by seqkit sana.
NEED_REUPLOAD  The gzip stream is incomplete or corrupted.
FAILED         Repair failed; inspect the log file.
```

The `input_for_skmer/` directory contains repaired files or links to valid originals and can be used as the workflow input directory.

## Outputs

All workflow outputs are written under `results/`.

### Cleaned Reads

```text
results/<sample>/clean/<sample>.R1.clean.fq.gz
results/<sample>/clean/<sample>.R2.clean.fq.gz
results/<sample>/clean/<sample>.fastp.html
results/<sample>/clean/<sample>.fastp.json
```

These are the `fastp` outputs. Use the HTML or JSON reports to inspect read quality, adapter trimming, filtering rates, and read length distributions.

### Reference-Filtered Reads

Produced only when `-ref` is supplied:

```text
results/<sample>/nDNA/un-conc-mate.1.fq.gz
results/<sample>/nDNA/un-conc-mate.2.fq.gz
results/<sample>/nDNA/<sample>.sam
```

The `un-conc-mate` files contain read pairs that did not align concordantly to the reference. These are retained for downstream analysis. The SAM file records alignments to the reference and can be used to inspect filtering behavior.

### Pair-Repaired Reads

```text
results/<sample>/nDNA/un-conc-mate.1.fixed.fq.gz
results/<sample>/nDNA/un-conc-mate.2.fixed.fq.gz
results/<sample>/nDNA/singletons.fq.gz
```

These files are generated after reference filtering. The `.fixed.fq.gz` files are synchronized paired reads. `singletons.fq.gz` contains reads whose mates were missing after filtering.

### Merged Reads

```text
results/<sample>/sample_merged.fq
results/<sample>/sample_unmerged1.fq
results/<sample>/sample_unmerged2.fq
results/<sample>/sample.fq
```

`sample_merged.fq` contains merged overlapping read pairs. The unmerged files contain pairs that could not be merged. `sample.fq` concatenates merged and unmerged reads for downstream statistics and normalization.

### Statistics And Normalized FASTQ

```text
results/stats/post_filter_summary.sorted.tsv
results/stats/head_base_cutoff.txt
results/stats/head_summary.sorted.tsv
results/<sample>/nDNAOK/<sample>.fq
```

`post_filter_summary.sorted.tsv` summarizes sequence depth after cleaning, filtering, merging, and concatenation.

`head_base_cutoff.txt` records the base-count cutoff selected by `-s`.

`head_summary.sorted.tsv` summarizes final normalized depth across samples.

`nDNAOK/<sample>.fq` is the final per-sample FASTQ used by Skmer, WASTER, and Mash.

### Skmer

```text
results/skmer/dimtrx_main.txt
results/skmer/dimtrx_main_cor_.txt
results/skmer/dimtrx_main_cor_OK.phy
results/skmer/tree.direct.tre
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
```

`dimtrx_main.txt` is the initial Skmer distance matrix. `dimtrx_main_cor_.txt` is the corrected distance matrix. `dimtrx_main_cor_OK.phy` is the PHYLIP-formatted matrix used by FastME.

`tree.direct.tre` is the direct tree. `tree.bootstrap.tre` is the bootstrap consensus tree. `tree.merged.tre` combines the direct tree with bootstrap support and is usually the main Skmer tree to inspect.

### WASTER

```text
results/waster/input.tsv
results/waster/waster.tree
```

`input.tsv` lists samples and final FASTQ paths. `waster.tree` is the WASTER tree inferred from the normalized reads.

### Mash

```text
results/mash/distances.tsv
results/mash/distances.phy
results/mash/distance_heatmap.svg
results/mash/tree.direct.tre
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
```

`distances.tsv` is the pairwise Mash distance matrix. `distance_heatmap.svg` visualizes these distances and is useful for detecting outliers or unexpected clustering.

`tree.direct.tre` is the direct Mash/FastME tree. `tree.bootstrap.tre` is the bootstrap consensus tree. `tree.merged.tre` combines the direct tree with bootstrap support.

## Recommended Files To Inspect

Quality control:

```text
results/<sample>/clean/<sample>.fastp.html
results/stats/post_filter_summary.sorted.tsv
results/stats/head_summary.sorted.tsv
results/mash/distance_heatmap.svg
```

Final trees:

```text
results/skmer/tree.merged.tre
results/mash/tree.merged.tre
results/waster/waster.tree
```

## Common Issues

### No Paired FASTQ Files Found

Make sure `-i` points to the FASTQ directory itself and that filenames use one of the supported R1/R2 naming styles.

### Missing External Tool

Run:

```bash
skmer-smk2 doctor
```

Then install the missing tool or add it to `PATH`.

### Existing Results Are Reused

Snakemake reuses completed outputs. For a clean rerun, use a new `--workdir` or remove old workflow outputs after confirming the paths.

### Old Workflow Cache Is Used

The workflow is materialized into:

```text
.skmer_smk2_workflow/
```

Remove this directory after upgrading if old commands still appear in `--printshellcmds`.
