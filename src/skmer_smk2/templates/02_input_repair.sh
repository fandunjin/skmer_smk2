#!/usr/bin/env bash
# Repair selected FASTQ samples after reviewing input_sample_report.tsv.

set -euo pipefail

# scheduler headers written by user

# source /path/to/conda.sh
# conda activate your_environment

FASTQ_DIR="${1:-/path/to/raw_fastq}"
OUT_DIR="${2:-input_qc}"
SAMPLES="${3:-sampleA sampleB}"

skmer-smk2 input-repair -i "${FASTQ_DIR}" -o "${OUT_DIR}" --samples "${SAMPLES}"

echo
echo "Review:"
echo "  ${OUT_DIR}/post_repair_sample_report.tsv"
echo
echo "Run skmer-smk2 with:"
echo "  skmer-smk2 run -i ${OUT_DIR}/input_for_skmer -ref /path/to/ref.fasta -s 75 -j 48"
