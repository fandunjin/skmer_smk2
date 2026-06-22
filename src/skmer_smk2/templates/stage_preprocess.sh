#!/bin/bash
# Run shared skmer-smk2 preprocessing only.
# Edit scheduler headers, environment activation, and paths for your cluster.

#JSUB -J skmer_prep
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_prep_log.%J
#JSUB -e skmer_prep_err.%J

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
  -prep \
  --fastp-threads 4 \
  --bowtie2-threads 2 \
  --repair-threads 2 \
  --bbmerge-threads 2 \
  --total-mem-mb 180000 \
  --bowtie2-mem-mb 6000 \
  --workdir "${RUN_DIR}" \
  --printshellcmds
