#!/bin/bash
#SBATCH --job-name=hisat2_for_braker
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hisat2_for_braker.%j.out
#SBATCH --error=hisat2_for_braker.%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
GENOME="/path/to/assembly_split.softmasked.fa"

READ1="/path/to/sample_R1.trimmed.fastq.gz"
READ2="/path/to/sample_R2.trimmed.fastq.gz"

OUTDIR="/path/to/output/hisat2_for_braker"
PREFIX="sample_name"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
INDEX_PREFIX="${OUTDIR}/genome_hisat2_index"

# -----------------------------
# Software requirements
# -----------------------------
# HISAT2
# samtools

# -----------------------------
# Checks
# -----------------------------
[[ -r "$GENOME" ]] || { echo "[ERR] Genome not found: $GENOME"; exit 1; }
[[ -r "$READ1" ]] || { echo "[ERR] READ1 not found: $READ1"; exit 1; }
[[ -r "$READ2" ]] || { echo "[ERR] READ2 not found: $READ2"; exit 1; }

command -v hisat2 >/dev/null || { echo "[ERR] hisat2 not in PATH"; exit 1; }
command -v hisat2-build >/dev/null || { echo "[ERR] hisat2-build not in PATH"; exit 1; }
command -v samtools >/dev/null || { echo "[ERR] samtools not in PATH"; exit 1; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -----------------------------
# Build HISAT2 index if missing
# -----------------------------
if [[ ! -f "${INDEX_PREFIX}.1.ht2" && ! -f "${INDEX_PREFIX}.1.ht2l" ]]; then
    echo "[INFO] Building HISAT2 index"
    hisat2-build -p "$THREADS" "$GENOME" "$INDEX_PREFIX"
else
    echo "[INFO] HISAT2 index already exists, skipping build"
fi

# -----------------------------
# Align RNA-seq reads
# -----------------------------
echo "[INFO] Running HISAT2"

hisat2 \
    -p "$THREADS" \
    --dta \
    -x "$INDEX_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    -S "${PREFIX}.sam" \
    2> "${PREFIX}.hisat2.stderr"

# -----------------------------
# Convert SAM to sorted BAM
# -----------------------------
echo "[INFO] Converting SAM to sorted BAM"

samtools view \
    -@ "$THREADS" \
    -b "${PREFIX}.sam" \
| samtools sort \
    -@ "$THREADS" \
    -o "${PREFIX}.sorted.bam"

# -----------------------------
# Index BAM
# -----------------------------
echo "[INFO] Indexing BAM"

samtools index -@ "$THREADS" "${PREFIX}.sorted.bam"

# -----------------------------
# Basic QC outputs
# -----------------------------
echo "[INFO] Writing flagstat and idxstats"

samtools flagstat -@ "$THREADS" "${PREFIX}.sorted.bam" > "${PREFIX}.flagstat.txt"
samtools idxstats "${PREFIX}.sorted.bam" > "${PREFIX}.idxstats.txt"

echo "[INFO] Saving BAM header"
samtools view -H "${PREFIX}.sorted.bam" > "${PREFIX}.header.txt"

echo "[INFO] Saving first 200 alignments for inspection"
samtools view "${PREFIX}.sorted.bam" | head -n 200 > "${PREFIX}.first200.txt" || true

rm -f "${PREFIX}.sam"

echo "[INFO] Done"
echo "[INFO] BAM: ${OUTDIR}/${PREFIX}.sorted.bam"
echo "[INFO] BAI: ${OUTDIR}/${PREFIX}.sorted.bam.bai"