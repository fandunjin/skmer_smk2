#!/bin/bash

set -euo pipefail

# Optional: activate the environment that contains skmer-smk2 and dependencies.
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate your_env

skmer-smk2 run \
  -i /path/to/fastq_dir \
  -ref /path/to/reference.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
