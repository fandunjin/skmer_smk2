#!/usr/bin/env bash
# Scan paired FASTQ files for skmer-smk2, repair usable record problems, and
# create a skmer-smk2-ready input directory.
#
# Usage:
#   bash scan_repair_fastq.sh /path/to/fastq_dir [output_dir]
#
# Outputs:
#   output_dir/fastq_repair_report.tsv
#   output_dir/sample_repair_report.tsv
#   output_dir/input_for_skmer/

set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    sed -n '1,12p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
fi

INPUT_DIR="${1:-$(pwd)}"
OUT_DIR="${2:-${INPUT_DIR}/repaired_fastq}"
SKMER_INPUT_DIR="${OUT_DIR}/input_for_skmer"
TMP_DIR="${OUT_DIR}/tmp"
REPORT="${OUT_DIR}/fastq_repair_report.tsv"
SAMPLE_REPORT="${OUT_DIR}/sample_repair_report.tsv"

mkdir -p "${OUT_DIR}" "${SKMER_INPUT_DIR}" "${TMP_DIR}"

need_tool() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: $1 was not found in PATH." >&2
        exit 1
    fi
}

need_tool gzip
need_tool seqkit
need_tool repair.sh
need_tool awk
need_tool sort

abs_path() {
    if command -v realpath >/dev/null 2>&1; then
        realpath "$1"
    else
        python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$1"
    fi
}

fastq_stem() {
    local name
    name="$(basename "$1")"
    name="${name%.gz}"
    name="${name%.fastq}"
    name="${name%.fq}"
    printf '%s\n' "${name}"
}

sample_from_stem() {
    local stem="$1"
    case "${stem}" in
        *_R1) printf '%s\n' "${stem%_R1}" ;;
        *_R2) printf '%s\n' "${stem%_R2}" ;;
        *.R1) printf '%s\n' "${stem%.R1}" ;;
        *.R2) printf '%s\n' "${stem%.R2}" ;;
        *-R1) printf '%s\n' "${stem%-R1}" ;;
        *-R2) printf '%s\n' "${stem%-R2}" ;;
        *_1) printf '%s\n' "${stem%_1}" ;;
        *_2) printf '%s\n' "${stem%_2}" ;;
        *.1) printf '%s\n' "${stem%.1}" ;;
        *.2) printf '%s\n' "${stem%.2}" ;;
        *-1) printf '%s\n' "${stem%-1}" ;;
        *-2) printf '%s\n' "${stem%-2}" ;;
        *) printf '%s\n' "" ;;
    esac
}

mate_from_stem() {
    local stem="$1"
    case "${stem}" in
        *_R1|*.R1|*-R1|*_1|*.1|*-1) printf 'R1\n' ;;
        *_R2|*.R2|*-R2|*_2|*.2|*-2) printf 'R2\n' ;;
        *) printf '\n' ;;
    esac
}

record_count() {
    seqkit stats -T "$1" 2>/dev/null | awk 'NR==2 {gsub(",", "", $4); print $4}'
}

base_count() {
    seqkit stats -T "$1" 2>/dev/null | awk 'NR==2 {gsub(",", "", $5); print $5}'
}

write_link_or_copy() {
    local src="$1"
    local dst="$2"
    rm -f "${dst}"
    if ln -s "$(abs_path "${src}")" "${dst}" 2>/dev/null; then
        return 0
    fi
    cp -f "${src}" "${dst}"
}

check_one_file() {
    local fq="$1"
    local base stem gzip_status seqkit_status records bases note
    base="$(basename "${fq}")"
    stem="$(fastq_stem "${fq}")"
    gzip_status="NA"
    seqkit_status="OK"
    records=""
    bases=""
    note=""

    if [[ "${fq}" == *.gz ]]; then
        if gzip -t "${fq}" 2>"${OUT_DIR}/${stem}.gzip.log"; then
            gzip_status="OK"
        else
            gzip_status="BAD"
            seqkit_status="SKIPPED"
            note="gzip stream failed; file is probably physically truncated and should be reuploaded"
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${base}" "${gzip_status}" "${seqkit_status}" "" "" "${note}" >> "${REPORT}"
            return 1
        fi
    fi

    if records="$(record_count "${fq}")" && bases="$(base_count "${fq}")" && [[ -n "${records}" ]]; then
        seqkit_status="OK"
        note="readable"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${base}" "${gzip_status}" "${seqkit_status}" "${records}" "${bases:-}" "${note}" >> "${REPORT}"
        return 0
    fi

    seqkit_status="BAD"
    note="seqkit stats failed"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${base}" "${gzip_status}" "${seqkit_status}" "" "" "${note}" >> "${REPORT}"
    return 2
}

sana_one_file() {
    local fq="$1"
    local out="$2"
    local stem
    stem="$(fastq_stem "${fq}")"
    seqkit sana "${fq}" -o "${out}" >"${OUT_DIR}/${stem}.seqkit_sana.log" 2>&1
}

printf 'file\tgzip_status\tseqkit_status\treads\tbases\tnote\n' > "${REPORT}"
printf 'sample\tr1\tr2\tr1_reads\tr2_reads\taction\toutput_r1\toutput_r2\tnote\n' > "${SAMPLE_REPORT}"

PAIR_TABLE="${TMP_DIR}/pairs.tsv"
: > "${PAIR_TABLE}"

find "${INPUT_DIR}" -maxdepth 1 -type f \( \
    -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fastq" \
\) | sort | while IFS= read -r fq; do
    stem="$(fastq_stem "${fq}")"
    sample="$(sample_from_stem "${stem}")"
    mate="$(mate_from_stem "${stem}")"
    if [[ -n "${sample}" && -n "${mate}" ]]; then
        printf '%s\t%s\t%s\n' "${sample}" "${mate}" "$(abs_path "${fq}")" >> "${PAIR_TABLE}"
    else
        printf '%s\tNA\tSKIPPED\t\t\tunsupported FASTQ mate naming style\n' "$(basename "${fq}")" >> "${REPORT}"
    fi
done

if [[ ! -s "${PAIR_TABLE}" ]]; then
    echo "ERROR: no paired FASTQ candidates found in ${INPUT_DIR}" >&2
    exit 1
fi

cut -f1 "${PAIR_TABLE}" | sort -u | while IFS= read -r sample; do
    r1="$(awk -F '\t' -v s="${sample}" '$1 == s && $2 == "R1" {print $3; exit}' "${PAIR_TABLE}")"
    r2="$(awk -F '\t' -v s="${sample}" '$1 == s && $2 == "R2" {print $3; exit}' "${PAIR_TABLE}")"

    if [[ -z "${r1}" || -z "${r2}" ]]; then
        printf '%s\t%s\t%s\t\t\tMISSING_PAIR\t\t\tmissing R1 or R2\n' "${sample}" "${r1:-}" "${r2:-}" >> "${SAMPLE_REPORT}"
        continue
    fi

    r1_ok=0
    r2_ok=0
    check_one_file "${r1}" || r1_ok=$?
    check_one_file "${r2}" || r2_ok=$?

    r1_reads="$(record_count "${r1}" || true)"
    r2_reads="$(record_count "${r2}" || true)"
    out_r1="${SKMER_INPUT_DIR}/$(basename "${r1}")"
    out_r2="${SKMER_INPUT_DIR}/$(basename "${r2}")"

    if [[ "${r1_ok}" -eq 0 && "${r2_ok}" -eq 0 && "${r1_reads}" == "${r2_reads}" && -n "${r1_reads}" ]]; then
        write_link_or_copy "${r1}" "${out_r1}"
        write_link_or_copy "${r2}" "${out_r2}"
        printf '%s\t%s\t%s\t%s\t%s\tUNCHANGED\t%s\t%s\tpassed gzip/readability/pair-count checks\n' \
            "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" "${r1_reads}" "${r2_reads}" "${out_r1}" "${out_r2}" >> "${SAMPLE_REPORT}"
        continue
    fi

    if [[ "${r1_ok}" -eq 1 || "${r2_ok}" -eq 1 ]]; then
        printf '%s\t%s\t%s\t%s\t%s\tNEED_REUPLOAD\t\t\tgzip failed before repair; reupload the damaged compressed file\n' \
            "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" "${r1_reads:-}" "${r2_reads:-}" >> "${SAMPLE_REPORT}"
        continue
    fi

    sana_r1="${OUT_DIR}/${sample}.R1.sana.fq.gz"
    sana_r2="${OUT_DIR}/${sample}.R2.sana.fq.gz"
    repaired_r1="${OUT_DIR}/${sample}.R1.repaired.fq.gz"
    repaired_r2="${OUT_DIR}/${sample}.R2.repaired.fq.gz"
    singletons="${OUT_DIR}/${sample}.singletons.fq.gz"
    repair_log="${OUT_DIR}/${sample}.repair.log"

    if ! sana_one_file "${r1}" "${sana_r1}" || ! sana_one_file "${r2}" "${sana_r2}"; then
        printf '%s\t%s\t%s\t%s\t%s\tSANA_FAILED\t\t\tseqkit sana failed; inspect *.seqkit_sana.log\n' \
            "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" "${r1_reads:-}" "${r2_reads:-}" >> "${SAMPLE_REPORT}"
        continue
    fi

    if repair.sh in1="${sana_r1}" in2="${sana_r2}" out1="${repaired_r1}" out2="${repaired_r2}" outs="${singletons}" repair=t overwrite=t >"${repair_log}" 2>&1; then
        repaired_r1_reads="$(record_count "${repaired_r1}" || true)"
        repaired_r2_reads="$(record_count "${repaired_r2}" || true)"
        if [[ -n "${repaired_r1_reads}" && "${repaired_r1_reads}" == "${repaired_r2_reads}" ]]; then
            write_link_or_copy "${repaired_r1}" "${out_r1}"
            write_link_or_copy "${repaired_r2}" "${out_r2}"
            printf '%s\t%s\t%s\t%s\t%s\tREPAIRED\t%s\t%s\trepaired_reads=%s; singletons=%s\n' \
                "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" "${r1_reads:-}" "${r2_reads:-}" \
                "${out_r1}" "${out_r2}" "${repaired_r1_reads}" "${singletons}" >> "${SAMPLE_REPORT}"
        else
            printf '%s\t%s\t%s\t%s\t%s\tREPAIR_FAILED\t\t\trepaired R1/R2 counts are empty or unequal\n' \
                "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" "${r1_reads:-}" "${r2_reads:-}" >> "${SAMPLE_REPORT}"
        fi
    else
        printf '%s\t%s\t%s\t%s\t%s\tREPAIR_FAILED\t\t\tBBMap repair.sh failed; inspect %s\n' \
            "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" "${r1_reads:-}" "${r2_reads:-}" "${repair_log}" >> "${SAMPLE_REPORT}"
    fi
done

echo
echo "FASTQ scan/repair finished."
echo "File report: ${REPORT}"
echo "Sample report: ${SAMPLE_REPORT}"
echo "Skmer-ready input directory: ${SKMER_INPUT_DIR}"
echo
echo "Sample action summary:"
awk -F '\t' 'NR > 1 {count[$6]++} END {for (k in count) print k, count[k]}' "${SAMPLE_REPORT}" | sort
echo
echo "Samples needing attention:"
awk -F '\t' 'NR == 1 || ($6 != "UNCHANGED" && $6 != "REPAIRED") {print}' "${SAMPLE_REPORT}"
