#!/usr/bin/env bash
#SBATCH --job-name=hisat2_edulis_144
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hisat2_ERR15315144_%j.out
#SBATCH --error=hisat2_ERR15315144_%j.err

set -euo pipefail

# Required software:
#   HISAT2
#   samtools

# -----------------------------
# Inputs
# -----------------------------

GENOME="/path/to/edulis_nuclear_only.softmasked.fa"

READ1="/path/to/ERR15315144_1_R1.trimmed.fastq.gz"
READ2="/path/to/ERR15315144_1_R2.trimmed.fastq.gz"

OUTDIR="/path/to/output/hisat2_edulis"
INDEX_PREFIX="${OUTDIR}/edulis_hisat2_index_ERR15315144"
PREFIX="ERR15315144_1"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# -----------------------------
# Build HISAT2 index
# -----------------------------

if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
    hisat2-build \
        -p "${THREADS}" \
        "${GENOME}" \
        "${INDEX_PREFIX}"
fi

# -----------------------------
# Align RNA-seq reads
# -----------------------------

hisat2 \
    -p "${THREADS}" \
    --dta \
    -x "${INDEX_PREFIX}" \
    -1 "${READ1}" \
    -2 "${READ2}" \
    -S "${PREFIX}.sam" \
    2> "${PREFIX}.hisat2.stderr"

# -----------------------------
# Convert to sorted BAM
# -----------------------------

samtools view \
    -@ "${THREADS}" \
    -bS "${PREFIX}.sam" \
| samtools sort \
    -@ "${THREADS}" \
    -o "${PREFIX}.sorted.bam"

samtools index \
    -@ "${THREADS}" \
    "${PREFIX}.sorted.bam"

samtools flagstat \
    -@ "${THREADS}" \
    "${PREFIX}.sorted.bam" \
    > "${PREFIX}.flagstat.txt"

samtools idxstats \
    "${PREFIX}.sorted.bam" \
    > "${PREFIX}.idxstats.txt"

rm -f "${PREFIX}.sam"

echo "[INFO] BAM: ${OUTDIR}/${PREFIX}.sorted.bam"
