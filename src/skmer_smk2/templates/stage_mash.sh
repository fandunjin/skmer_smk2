#!/bin/bash
# Run only the Mash branch after preprocessing.
# Mash contains many small single-process jobs, so a high -j value is useful.

#JSUB -J mash_only
#JSUB -n 48
#JSUB -q normal
#JSUB -o mash_only_log.%J
#JSUB -e mash_only_err.%J

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
  -j 48 \
  -mash \
  --total-mem-mb 80000 \
  --workdir "${RUN_DIR}" \
  --printshellcmds
