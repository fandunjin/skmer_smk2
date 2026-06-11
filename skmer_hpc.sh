#!/bin/bash
#JSUB -J skmer
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_log.%J
#JSUB -e skmer_err.%J
#JSUB -cwd skmer_%J

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
EXCLUDE_NORMALIZED="${EXCLUDE_SAMPLES//,/ }"

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

test -f ./snakefile
test -d ./scripts
test -f "${REF}"

for tool in snakemake python fastp bowtie2 repair.sh bbmerge.sh skmer fastme raxmlHPC waster mash; do
    command -v "${tool}" >/dev/null 2>&1 || {
        echo "ERROR: ${tool} was not found in the active conda environment." >&2
        exit 1
    }
done

echo "Checking FASTQ gzip integrity..."
BAD_FASTQ=0
for f in ./*_1.fq.gz ./*_2.fq.gz; do
    test -e "${f}" || continue
    base=$(basename "${f}")
    sample="${base%_1.fq.gz}"
    sample="${sample%_2.fq.gz}"
    if [[ " ${EXCLUDE_NORMALIZED} " == *" ${sample} "* ]]; then
        echo "Skipping excluded sample FASTQ: ${f}"
        continue
    fi
    gzip -t "${f}" || {
        echo "BAD FASTQ gzip: ${f}" >&2
        BAD_FASTQ=1
    }
done
if [[ "${BAD_FASTQ}" -ne 0 ]]; then
    echo "ERROR: at least one non-excluded FASTQ gzip file is incomplete or corrupted." >&2
    echo "Re-upload the BAD FASTQ gzip files, or add the sample name to EXCLUDE_SAMPLES." >&2
    exit 1
fi
echo "FASTQ gzip integrity check passed."

snakemake \
    -s ./snakefile \
    --cores "${THREADS}" \
    --rerun-incomplete \
    --latency-wait 120 \
    --printshellcmds \
    --config \
        input_dir="${INPUT_DIR}" \
        ref="${REF}" \
        sample_percentile="${SAMPLE_PERCENTILE}" \
        rep_n="${BOOTSTRAPS}" \
        exclude_samples="${EXCLUDE_SAMPLES}"

echo
echo "Job finished at: $(date)"
