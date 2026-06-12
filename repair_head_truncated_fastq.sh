#!/bin/bash
#JSUB -J repair_fastq
#JSUB -n 8
#JSUB -q normal
#JSUB -cwd .
#JSUB -o repair_fastq.%J.out
#JSUB -e repair_fastq.%J.err

set -u
set -o pipefail

THREADS="${THREADS:-8}"
OUTDIR="${OUTDIR:-repaired_for_skmer}"
TMPDIR="${TMPDIR:-repair_tmp}"
REPORT="${REPORT:-fastq_repair_report.tsv}"
BADDIR="${BADDIR:-failed_repair_outputs}"

# Some conda activate.d scripts reference unset variables, which breaks under
# `set -u`. Temporarily relax nounset only while activating the environment.
_nounset_was_on=0
[[ $- == *u* ]] && _nounset_was_on=1
set +u
source /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh
conda activate 01bio
if (( _nounset_was_on )); then
  set -u
fi
unset _nounset_was_on

samples=(
  "H_miyoshiana_SAMC1020838"
  "H_serrata_SAMC1020837"
  "H_asiatica_SAMC1020836"
  "H_lucidula_SAMC1020839"
)

mkdir -p "$OUTDIR" "$TMPDIR" "$BADDIR"
printf "sample\tr1_input\tr2_input\tr1_output\tr2_output\tstatus\tnote\n" > "$REPORT"

need_tool() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: required tool not found in PATH: $1" >&2
    exit 127
  }
}

need_tool seqkit
need_tool gzip

fastq_reads() {
  seqkit stats -T "$1" | awk 'NR==2 {gsub(/,/, "", $4); print $4}'
}

quarantine_outputs() {
  local sample="$1"
  shift
  mkdir -p "${BADDIR}/${sample}"
  for path in "$@"; do
    if [[ -e "$path" ]]; then
      mv -f "$path" "${BADDIR}/${sample}/"
    fi
  done
}

for sample in "${samples[@]}"; do
  r1="${sample}_1.fq.gz"
  r2="${sample}_2.fq.gz"
  s1="${TMPDIR}/${sample}_1.sana.fq.gz"
  s2="${TMPDIR}/${sample}_2.sana.fq.gz"
  pairdir="${TMPDIR}/${sample}.pair.$$"
  o1="${OUTDIR}/${sample}_1.fq.gz"
  o2="${OUTDIR}/${sample}_2.fq.gz"
  log="${TMPDIR}/${sample}.seqkit.log"

  echo "============================================================"
  echo "Processing ${sample}"

  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\tMISSING_INPUT\tinput file not found or empty\n" "$sample" "$r1" "$r2" "$o1" "$o2" >> "$REPORT"
    echo "WARNING: missing input for ${sample}" >&2
    continue
  fi

  {
    echo "[sana R1] $r1"
    seqkit sana "$r1" -o "$s1"
    echo "[sana R2] $r2"
    seqkit sana "$r2" -o "$s2"
    echo "[pair] $sample"
    seqkit pair -1 "$s1" -2 "$s2" -O "$pairdir" -j "$THREADS"
  } > "$log" 2>&1

  mv -f "${pairdir}/$(basename "$s1")" "$o1" 2>/dev/null || true
  mv -f "${pairdir}/$(basename "$s2")" "$o2" 2>/dev/null || true

  if [[ ! -s "$o1" || ! -s "$o2" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\tFAILED\tseqkit output missing; see %s\n" "$sample" "$r1" "$r2" "$o1" "$o2" "$log" >> "$REPORT"
    echo "ERROR: output missing for ${sample}; see ${log}" >&2
    continue
  fi

  if gzip -t "$o1" >/dev/null 2>&1 && gzip -t "$o2" >/dev/null 2>&1; then
    reads1=$(fastq_reads "$o1")
    reads2=$(fastq_reads "$o2")
    if [[ ! "$reads1" =~ ^[0-9]+$ || ! "$reads2" =~ ^[0-9]+$ ]]; then
      printf "%s\t%s\t%s\t%s\t%s\tBAD_OUTPUT\tcould not read FASTQ counts; see %s\n" "$sample" "$r1" "$r2" "$o1" "$o2" "$log" >> "$REPORT"
      echo "ERROR: could not read FASTQ counts for ${sample}; see ${log}" >&2
      quarantine_outputs "$sample" "$o1" "$o2"
      continue
    fi
    if (( reads1 == 0 || reads2 == 0 )); then
      printf "%s\t%s\t%s\t%s\t%s\tEMPTY_OUTPUT\tpaired output has zero reads and was moved to %s/%s; reads R1=%s R2=%s; see %s\n" "$sample" "$r1" "$r2" "$o1" "$o2" "$BADDIR" "$sample" "$reads1" "$reads2" "$log" >> "$REPORT"
      echo "ERROR: paired output has zero reads for ${sample}; see ${log}" >&2
      quarantine_outputs "$sample" "$o1" "$o2"
      continue
    fi
    if (( reads1 != reads2 )); then
      printf "%s\t%s\t%s\t%s\t%s\tUNPAIRED_OUTPUT\tpaired output read counts differ; reads R1=%s R2=%s; see %s\n" "$sample" "$r1" "$r2" "$o1" "$o2" "$reads1" "$reads2" "$log" >> "$REPORT"
      echo "ERROR: paired output read counts differ for ${sample}; see ${log}" >&2
      quarantine_outputs "$sample" "$o1" "$o2"
      continue
    fi
    printf "%s\t%s\t%s\t%s\t%s\tREPAIRED\tpaired output gzip OK; reads R1=%s R2=%s\n" "$sample" "$r1" "$r2" "$o1" "$o2" "$reads1" "$reads2" >> "$REPORT"
    echo "OK: ${sample}"
    rm -f "$s1" "$s2"
    rmdir "$pairdir" 2>/dev/null || true
  else
    printf "%s\t%s\t%s\t%s\t%s\tBAD_OUTPUT\tgzip -t failed after repair; see %s\n" "$sample" "$r1" "$r2" "$o1" "$o2" "$log" >> "$REPORT"
    echo "ERROR: gzip test failed after repair for ${sample}" >&2
    quarantine_outputs "$sample" "$o1" "$o2"
  fi
done

echo "============================================================"
echo "Done."
echo "Repaired FASTQ directory: $OUTDIR"
echo "Report: $REPORT"
echo
echo "Next step example:"
echo "  skmer-smk2 run -i $OUTDIR -ref /path/to/refDNA.fasta -s 75 -j 48 --printshellcmds"
