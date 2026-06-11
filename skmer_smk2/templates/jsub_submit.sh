#!/bin/bash
#JSUB -J skmer_smk2
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_smk2.%J.out
#JSUB -e skmer_smk2.%J.err
#JSUB -cwd skmer_smk2_%J

set -euo pipefail

WORKDIR="${WORKDIR:-/path/to/workdir}"
INPUT_DIR="${INPUT_DIR:-${WORKDIR}}"
REF="${REF:-${WORKDIR}/ref/refDNA.fasta}"
THREADS="${THREADS:-48}"
SAMPLE_PERCENTILE="${SAMPLE_PERCENTILE:-75}"
BOOTSTRAPS="${BOOTSTRAPS:-100}"
EXCLUDE_SAMPLES="${EXCLUDE_SAMPLES:-}"
CONDA_PROFILE="${CONDA_PROFILE:-}"
CONDA_ENV="${CONDA_ENV:-}"

if [[ -n "${CONDA_PROFILE}" ]]; then
    source "${CONDA_PROFILE}"
fi

if [[ -n "${CONDA_ENV}" ]]; then
    set +u
    conda activate "${CONDA_ENV}"
    set -u
fi

cd "${WORKDIR}"

echo "Job started at: $(date)"
echo "Working directory: ${WORKDIR}"
echo "Input FASTQ directory: ${INPUT_DIR}"
echo "Reference genome: ${REF}"
echo "Threads: ${THREADS}"
echo "Sample percentile: ${SAMPLE_PERCENTILE}"
echo "Bootstrap replicates: ${BOOTSTRAPS}"
echo "Excluded samples: ${EXCLUDE_SAMPLES:-none}"
echo

for tool in snakemake python fastp bowtie2 repair.sh bbmerge.sh skmer fastme raxmlHPC waster mash skmer-smk2; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: ${tool} was not found in PATH." >&2
        exit 1
    }
done

skmer-smk2 run \
    --workdir "${WORKDIR}" \
    -i "${INPUT_DIR}" \
    -ref "${REF}" \
    -s "${SAMPLE_PERCENTILE}" \
    -j "${THREADS}" \
    -b "${BOOTSTRAPS}" \
    --exclude-samples "${EXCLUDE_SAMPLES}" \
    --printshellcmds

echo
echo "Job finished at: $(date)"
