#!/bin/bash
# Run only the Skmer branch after preprocessing.
# Existing preprocessing outputs are reused automatically.

#JSUB -J skmer_only
#JSUB -n 16
#JSUB -q normal
#JSUB -o skmer_only_log.%J
#JSUB -e skmer_only_err.%J

set -euo pipefail

# source /path/to/conda/etc/profile.d/conda.sh
# conda activate your_env

FASTQ_DIR="/path/to/fastq_dir"
REF_FASTA="/path/to/ref.fasta"  # Leave empty when no reference filtering is needed.
RUN_DIR="/path/to/run_directory"

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

REF_ARGS=()
if [ -n "${REF_FASTA}" ]; then
  REF_ARGS=(-ref "${REF_FASTA}")
fi

skmer-smk2 run \
  -i "${FASTQ_DIR}" \
  "${REF_ARGS[@]}" \
  -s 75 \
  --candidate-percentiles 40,50,60,75,90 \
  -j 16 \
  -skmer \
  --skmer-sketch-size 50000 \
  --skmer-threads 8 \
  --skmer-mem-mb 180000 \
  --total-mem-mb 200000 \
  --workdir "${RUN_DIR}" \
  --printshellcmds
