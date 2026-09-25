#!/bin/bash
#SBATCH --job-name=eggnog_edulis
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=eggnog_edulis_%j.out
#SBATCH --error=eggnog_edulis_%j.err

set -euo pipefail

# Required software:
#   eggNOG-mapper
#   DIAMOND

# -----------------------------
# Inputs
# -----------------------------

INPUT_FASTA="/path/to/edulis_braker3/braker.longest_isoform.gffread.aa"

OUTDIR="/path/to/output/eggnog_edulis"
PREFIX="edulis"

EGGNOG_DATA_DIR="/path/to/eggnog_db"

THREADS="${SLURM_CPUS_PER_TASK:-24}"

mkdir -p "${OUTDIR}"

# Use scheduler-provided temporary storage when available.
TMP="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}"

# -----------------------------
# Run eggNOG-mapper
# -----------------------------

emapper.py \
    -i "${INPUT_FASTA}" \
    --itype proteins \
    -m diamond \
    --cpu "${THREADS}" \
    --data_dir "${EGGNOG_DATA_DIR}" \
    --output "${PREFIX}" \
    --output_dir "${OUTDIR}" \
    --temp_dir "${TMP}" \
    --override \
    --no_file_comments \
    --report_orthologs

echo "[INFO] Output: ${OUTDIR}/${PREFIX}.emapper.annotations"
