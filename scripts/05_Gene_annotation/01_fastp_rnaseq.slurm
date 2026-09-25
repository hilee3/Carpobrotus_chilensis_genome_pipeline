#!/bin/bash
#SBATCH --job-name=fastp_rnaseq
#SBATCH --output=fastp_%j.out
#SBATCH --error=fastp_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -euo pipefail

# -----------------------------
# Input
# -----------------------------
# C. chilensis RNA-seq data
# NCBI SRA accession: SRR38839563
# Edit these paths before running.

R1="/path/to/raw/sample_R1.fastq.gz"
R2="/path/to/raw/sample_R2.fastq.gz"
SAMPLE="sample_name"

# -----------------------------
# Output
# -----------------------------
OUTDIR="/path/to/output/trimmed_fastp"
REPORTDIR="${OUTDIR}/reports"

mkdir -p "$OUTDIR" "$REPORTDIR"

# -----------------------------
# Software requirements
# -----------------------------
# fastp

command -v fastp >/dev/null || { echo "[ERR] fastp not in PATH"; exit 1; }

[[ -r "$R1" ]] || { echo "[ERR] R1 not found: $R1"; exit 1; }
[[ -r "$R2" ]] || { echo "[ERR] R2 not found: $R2"; exit 1; }

# -----------------------------
# Adapter and quality trimming
# -----------------------------
echo "[INFO] Running fastp for ${SAMPLE}"

fastp \
    -i "$R1" \
    -I "$R2" \
    -o "${OUTDIR}/${SAMPLE}_R1.trimmed.fastq.gz" \
    -O "${OUTDIR}/${SAMPLE}_R2.trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --length_required 35 \
    --thread "${SLURM_CPUS_PER_TASK:-8}" \
    --overrepresentation_analysis \
    --html "${REPORTDIR}/${SAMPLE}.fastp.html" \
    --json "${REPORTDIR}/${SAMPLE}.fastp.json"

echo "[INFO] fastp preprocessing finished for ${SAMPLE}"