#!/bin/bash
#SBATCH --job-name=eggnog_mapper
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=eggnog_%j.out
#SBATCH --error=eggnog_%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.

# Longest-isoform protein FASTA generated from the BRAKER3 annotation.
# Transcript g103.t2 was removed prior to longest-isoform selection
# because of an internal-stop issue.
INPUT_FASTA="/path/to/braker.g103t2_removed.longest_isoform.gffread.aa"

# Output prefix. The parent directory will be created if it does not exist.
OUT_PREFIX="/path/to/output/eggnog/output"

# eggNOG database directory.
export EGGNOG_DATA_DIR="/path/to/eggnog_db"

# Optional temp/scratch directory for large intermediates.
if [[ -z "${TMPDIR:-}" ]]; then
    export TMPDIR="${SLURM_TMPDIR:-/tmp}"
fi

# -----------------------------
# Software requirements
# -----------------------------
# eggNOG-mapper must be available in PATH.
# This script was run using DIAMOND mode with protein input.

command -v emapper.py >/dev/null || {
    echo "[ERR] emapper.py not found in PATH"
    exit 1
}

[[ -r "$INPUT_FASTA" ]] || {
    echo "[ERR] Input protein FASTA not found: $INPUT_FASTA"
    exit 1
}

[[ -d "$EGGNOG_DATA_DIR" ]] || {
    echo "[ERR] EGGNOG_DATA_DIR not found: $EGGNOG_DATA_DIR"
    exit 1
}

OUT_DIR="$(dirname "$OUT_PREFIX")"
OUT_NAME="$(basename "$OUT_PREFIX")"
mkdir -p "$OUT_DIR"

# -----------------------------
# Run eggNOG-mapper
# -----------------------------
echo "[INFO] Starting eggNOG-mapper: $(date)"
echo "[INFO] Input FASTA: $INPUT_FASTA"
echo "[INFO] Data dir: $EGGNOG_DATA_DIR"
echo "[INFO] Output dir: $OUT_DIR"
echo "[INFO] Output name: $OUT_NAME"
echo "[INFO] TMPDIR: $TMPDIR"

emapper.py \
    -i "$INPUT_FASTA" \
    --itype proteins \
    -m diamond \
    --cpu "${SLURM_CPUS_PER_TASK:-24}" \
    --data_dir "$EGGNOG_DATA_DIR" \
    --output "$OUT_NAME" \
    --output_dir "$OUT_DIR" \
    --temp_dir "$TMPDIR" \
    --override \
    --no_file_comments \
    --report_orthologs

echo "[INFO] eggNOG-mapper finished: $(date)"

# -----------------------------
# Quick summary
# -----------------------------
ANNOT="${OUT_PREFIX}.emapper.annotations"

if [[ -f "$ANNOT" ]]; then
    echo "[INFO] Summary of annotated proteins:"
    awk 'BEGIN{FS=OFS="\t"} $1!~/^#/ {n++; if($6!="-" || $8!="-" || $10!="-") ann++} END{print "Total_rows:", n, "Annotated_rows:", ann+0}' "$ANNOT"
else
    echo "[WARN] Annotation file not found: $ANNOT"
fi
