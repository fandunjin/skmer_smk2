#!/bin/bash
# Scan FASTQ files, repair malformed records with seqkit sana, and write a report.
#
# Usage:
#   bash scan_repair_fastq.sh
#   bash scan_repair_fastq.sh /path/to/fastq_dir
#
# Outputs:
#   repaired_fastq/fastq_repair_report.tsv
#   repaired_fastq/*.repaired.fq.gz for files that seqkit sana repaired
#   repaired_fastq/input_for_skmer/ with repaired files or symlinks to good originals

set -euo pipefail

INPUT_DIR="${1:-$(pwd)}"
OUT_DIR="${INPUT_DIR}/repaired_fastq"
SKMER_INPUT_DIR="${OUT_DIR}/input_for_skmer"
REPORT="${OUT_DIR}/fastq_repair_report.tsv"

mkdir -p "${OUT_DIR}" "${SKMER_INPUT_DIR}"

if ! command -v seqkit >/dev/null 2>&1; then
    echo "ERROR: seqkit was not found in PATH. Activate the environment containing seqkit first." >&2
    exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip was not found in PATH." >&2
    exit 1
fi

printf "file\tgzip_status\tseqkit_status\taction\toutput\tnote\n" > "${REPORT}"

find "${INPUT_DIR}" -maxdepth 1 -type f \( \
    -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fastq" \
\) | sort | while IFS= read -r fq; do
    base=$(basename "${fq}")
    stem="${base}"
    stem="${stem%.gz}"
    stem="${stem%.fastq}"
    stem="${stem%.fq}"

    repaired="${OUT_DIR}/${stem}.repaired.fq.gz"
    sana_log="${OUT_DIR}/${stem}.seqkit_sana.log"
    skmer_link="${SKMER_INPUT_DIR}/${base}"

    gzip_status="NA"
    if [[ "${fq}" == *.gz ]]; then
        if gzip -t "${fq}" 2>"${OUT_DIR}/${stem}.gzip.log"; then
            gzip_status="OK"
        else
            gzip_status="BAD"
            printf "%s\t%s\tSKIPPED\tNEED_REUPLOAD\t\tgzip stream is incomplete/corrupted; seqkit cannot reliably repair truncated compressed files\n" \
                "${base}" "${gzip_status}" >> "${REPORT}"
            continue
        fi
    fi

    if seqkit sana "${fq}" -o "${repaired}" >"${sana_log}" 2>&1; then
        if seqkit stats -T "${repaired}" >/dev/null 2>&1; then
            seqkit_status="OK"
        else
            seqkit_status="BAD_AFTER_SANA"
            printf "%s\t%s\t%s\tFAILED\t%s\tseqkit sana finished but repaired file failed seqkit stats\n" \
                "${base}" "${gzip_status}" "${seqkit_status}" "${repaired}" >> "${REPORT}"
            continue
        fi
    else
        seqkit_status="SANA_FAILED"
        printf "%s\t%s\t%s\tFAILED\t%s\tsee %s\n" \
            "${base}" "${gzip_status}" "${seqkit_status}" "${repaired}" "${sana_log}" >> "${REPORT}"
        continue
    fi

    original_records=$(seqkit stats -T "${fq}" 2>/dev/null | awk 'NR==2 {gsub(",", "", $4); print $4}')
    repaired_records=$(seqkit stats -T "${repaired}" 2>/dev/null | awk 'NR==2 {gsub(",", "", $4); print $4}')

    if [[ "${original_records:-}" == "${repaired_records:-}" ]]; then
        rm -f "${repaired}"
        ln -sf "$(realpath "${fq}")" "${skmer_link}"
        printf "%s\t%s\t%s\tUNCHANGED\t%s\trecords=%s\n" \
            "${base}" "${gzip_status}" "${seqkit_status}" "${skmer_link}" "${original_records:-unknown}" >> "${REPORT}"
    else
        ln -sf "$(realpath "${repaired}")" "${skmer_link}"
        printf "%s\t%s\t%s\tREPAIRED\t%s\toriginal_records=%s; repaired_records=%s\n" \
            "${base}" "${gzip_status}" "${seqkit_status}" "${repaired}" "${original_records:-unknown}" "${repaired_records:-unknown}" >> "${REPORT}"
    fi
done

echo
echo "FASTQ scan/repair finished."
echo "Report: ${REPORT}"
echo "Skmer-ready input directory: ${SKMER_INPUT_DIR}"
echo
echo "Summary:"
awk -F '\t' 'NR > 1 {count[$4]++} END {for (k in count) print k, count[k]}' "${REPORT}" | sort
echo
echo "Problematic files:"
awk -F '\t' 'NR == 1 || $4 != "UNCHANGED" {print}' "${REPORT}"
