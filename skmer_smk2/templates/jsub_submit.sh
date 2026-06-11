#!/bin/bash
#JSUB -J skmer_smk2
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_smk2.%J.out
#JSUB -e skmer_smk2.%J.err

set -euo pipefail

# Edit these two lines for your own environment.
source /path/to/conda/etc/profile.d/conda.sh
conda activate your_env

# The workflow command is the same on local machines and HPC compute nodes.
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/refDNA.fasta -s 75 -j 48 --printshellcmds
