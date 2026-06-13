#!/usr/bin/env bash
# Check FASTQ inputs before running skmer-smk2.

set -euo pipefail

# scheduler headers written by user

# source /path/to/conda.sh
# conda activate your_environment

FASTQ_DIR="${1:-/path/to/raw_fastq}"
OUT_DIR="${2:-input_qc}"

skmer-smk2 input-check -i "${FASTQ_DIR}" -o "${OUT_DIR}"

echo
echo "Review:"
echo "  ${OUT_DIR}/input_sample_report.tsv"
