# skmer_snakemake_project

Snakemake workflow for paired-end FASTQ processing, plastid genome read removal, base-aware subsampling, and phylogenetic tree inference with Skmer, WASTER, and Mash.

This project is an updated workflow derived from the style and example-data layout of [`fandunjin/skmer_smk`](https://github.com/fandunjin/skmer_smk). The bundled `raw_data` demo reference and FASTQ files are copied from that repository.

## Workflow

The workflow performs:

1. Paired-end FASTQ discovery from `_1/_2` or `_R1/_R2` naming.
2. Read quality filtering with `fastp`.
3. Optional plastid genome removal with `bowtie2` when a reference is supplied.
4. Pair repair with `repair.sh`.
5. Read merging with `bbmerge.sh`.
6. FASTQ statistics after filtering and after final subsampling.
7. Base-aware subsampling using the sorted sample-depth percentile.
8. Skmer direct tree, bootstrap tree, and merged tree.
9. WASTER tree.
10. Mash direct tree, bootstrap tree, merged tree, and distance heatmap.

## Project Structure

```text
.
├── snakefile
├── run_skmer.py
├── skmer_hpc.sh
├── scan_repair_fastq.sh
├── raw_data/
│   ├── ref.fna
│   └── raw_data/
│       ├── sample1.R1.fq.gz
│       ├── sample1.R2.fq.gz
│       └── ...
├── scripts/
│   ├── distance_to_phylip.py
│   ├── fastq_stats_and_sample.py
│   ├── mash_dist_to_phylip.py
│   ├── merge_consensus.py
│   ├── plot_mash_heatmap.py
│   ├── summarize_fastq_stats.py
│   └── write_sample_fastq_list.py
└── README.md
```

Large FASTQ files, reference files, and workflow outputs are intentionally excluded from Git.

## Input Data

Use paired-end FASTQ files in one directory:

```text
SampleA_1.fq.gz
SampleA_2.fq.gz
SampleB_R1.fq.gz
SampleB_R2.fq.gz
SampleC.R1.fq.gz
SampleC.R2.fq.gz
```

For plastid removal, provide a reference FASTA:

```text
ref/refDNA.fasta
```

Bundled demo data:

```text
raw_data/ref.fna
raw_data/raw_data/sample*.R1.fq.gz
raw_data/raw_data/sample*.R2.fq.gz
```

## Dependencies

The active environment should include:

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

## Run

Without plastid filtering:

```bash
python run_skmer.py -i /path/to/fastq_dir -s 75 -j 48
```

With plastid filtering:

```bash
python run_skmer.py -i /path/to/fastq_dir -ref /path/to/refDNA.fasta -s 75 -j 48
```

Demo command:

```bash
python run_skmer.py -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 4 --dry-run
```

The bundled demo data has been checked with a local Snakemake dry-run using `rep_n=2`. Full execution requires the bioinformatics tools listed above, preferably in the HPC conda environment.

The value of `-s` is the sorted sample-depth percentile used to choose the base-count cutoff. For example, `-s 75` sorts samples by total bases from large to small, takes the base count at the 75% position, and uses that value for final FASTQ truncation.

## HPC Example

Edit paths in `skmer_hpc.sh`, then submit:

```bash
cd /hpcfile/users/92024286/Huperzia
jsub < skmer_hpc.sh
```

Important variables in `skmer_hpc.sh`:

```bash
WORKDIR=/hpcfile/users/92024286/Huperzia
INPUT_DIR="${WORKDIR}"
REF="${WORKDIR}/ref/refDNA.fasta"
THREADS=48
SAMPLE_PERCENTILE=75
BOOTSTRAPS=100
EXCLUDE_SAMPLES="H_serrata_SAMC1020837"
```

## FASTQ Scan And Repair

Before running the full workflow, scan FASTQ integrity and repair malformed FASTQ records:

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
