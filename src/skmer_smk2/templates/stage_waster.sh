#!/bin/bash
# Run only the WASTER branch after preprocessing.
# WASTER is usually one large job, so request few cores but enough memory.

#JSUB -J waster_only
#JSUB -n 1
#JSUB -q normal
#JSUB -o waster_only_log.%J
#JSUB -e waster_only_err.%J

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
  -j 1 \
  -waster \
  --waster-threads 1 \
  --waster-mem-mb 220000 \
  --total-mem-mb 240000 \
  --workdir "${RUN_DIR}" \
  --printshellcmds
