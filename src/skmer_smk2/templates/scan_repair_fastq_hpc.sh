#!/usr/bin/env bash
# HPC wrapper for the skmer-smk2 FASTQ scan/repair helper.
#
# Edit only the scheduler header and environment activation section for your
# cluster, then submit this script with your scheduler command.

# scheduler headers written by user

set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    sed -n '1,9p' "$0" | sed 's/^# \{0,1\}//'
    echo
    echo "Usage: bash scan_repair_fastq_hpc.sh /path/to/fastq_dir [output_dir]"
    exit 0
fi

# Example environment activation. Edit for your system.
# source /path/to/conda.sh
# conda activate your_environment

FASTQ_DIR="${1:-$(pwd)}"
OUT_DIR="${2:-${FASTQ_DIR}/repaired_fastq}"

if command -v skmer-smk2 >/dev/null 2>&1; then
    skmer-smk2 repair-fastq -i "${FASTQ_DIR}" --workdir "${OUT_DIR}" --copy-only
    bash "${OUT_DIR}/scan_repair_fastq.sh" "${FASTQ_DIR}" "${OUT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    bash "${SCRIPT_DIR}/scan_repair_fastq.sh" "${FASTQ_DIR}" "${OUT_DIR}"
fi
