#!/bin/bash
#JSUB -J replace_repaired_fastq
#JSUB -n 1
#JSUB -q normal
#JSUB -cwd .
#JSUB -o replace_repaired_fastq.%J.out
#JSUB -e replace_repaired_fastq.%J.err

set -euo pipefail

OUTDIR="${OUTDIR:-repaired_for_skmer}"
BACKUP_DIR="${BACKUP_DIR:-original_fastq_backup_$(date +%Y%m%d_%H%M%S)}"

samples=(
  "H_asiatica_SAMC1020836"
  "H_lucidula_SAMC1020839"
)

_nounset_was_on=0
[[ $- == *u* ]] && _nounset_was_on=1
set +u
if [[ -f /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh ]]; then
  source /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh
  conda activate 01bio
fi
if (( _nounset_was_on )); then
  set -u
fi
unset _nounset_was_on

need_tool() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: required tool not found in PATH: $1" >&2
    exit 127
  }
}

fastq_reads() {
  seqkit stats -T "$1" | awk 'NR==2 {gsub(/,/, "", $4); print $4}'
}

need_tool gzip
need_tool seqkit

echo "Backup directory: ${BACKUP_DIR}"
mkdir -p "$BACKUP_DIR"

for sample in "${samples[@]}"; do
  src1="${sample}_1.fq.gz"
  src2="${sample}_2.fq.gz"
  repaired1="${OUTDIR}/${sample}_1.fq.gz"
  repaired2="${OUTDIR}/${sample}_2.fq.gz"

  echo "============================================================"
  echo "Checking ${sample}"

  for path in "$src1" "$src2" "$repaired1" "$repaired2"; do
    if [[ ! -s "$path" ]]; then
      echo "ERROR: missing or empty file: $path" >&2
      exit 1
    fi
  done

  gzip -t "$repaired1"
  gzip -t "$repaired2"

  reads1=$(fastq_reads "$repaired1")
  reads2=$(fastq_reads "$repaired2")
  if [[ ! "$reads1" =~ ^[0-9]+$ || ! "$reads2" =~ ^[0-9]+$ ]]; then
    echo "ERROR: could not read FASTQ counts for ${sample}" >&2
    exit 1
  fi
  if (( reads1 == 0 || reads2 == 0 || reads1 != reads2 )); then
    echo "ERROR: repaired files are not a valid paired set for ${sample}: R1=${reads1}, R2=${reads2}" >&2
    exit 1
  fi

  mv -n "$src1" "${BACKUP_DIR}/"
  mv -n "$src2" "${BACKUP_DIR}/"
  cp -p "$repaired1" "$src1"
  cp -p "$repaired2" "$src2"

  echo "Replaced ${sample}: reads R1=${reads1} R2=${reads2}"
done

echo "============================================================"
echo "Done."
echo "Original FASTQ backup: ${BACKUP_DIR}"
echo "Repaired files copied from: ${OUTDIR}"
