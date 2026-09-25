#!/bin/bash
#SBATCH --job-name=hifi_mapping
#SBATCH --output=hifi_mapping_%j.out
#SBATCH --error=hifi_mapping_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G

set -euo pipefail

# Required software:
#   minimap2
#   samtools

# -----------------------------
# Inputs
# -----------------------------

REF="/path/to/assembly_split.fa"
READS="/path/to/pacbio_hifi_reads.fasta.gz"

OUTDIR="/path/to/output"
PREFIX="hifi_cov_split"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"


# -----------------------------
# Index reference
# -----------------------------

if [[ ! -s "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi


# -----------------------------
# Map HiFi reads to assembly
# -----------------------------

minimap2 \
    -t "${THREADS}" \
    -ax map-hifi \
    "${REF}" \
    "${READS}" \
| samtools view \
    -@ "${THREADS}" \
    -b \
    - \
| samtools sort \
    -@ "${THREADS}" \
    -o "${PREFIX}_sorted.bam" \
    -


# -----------------------------
# Index BAM
# -----------------------------

# CSI indexing was used for the coverage BAM.
samtools index \
    -c \
    -@ "${THREADS}" \
    "${PREFIX}_sorted.bam"


# -----------------------------
# Basic alignment summaries
# -----------------------------

samtools quickcheck \
    "${PREFIX}_sorted.bam"

samtools idxstats \
    "${PREFIX}_sorted.bam" \
    > "${PREFIX}.idxstats.txt"

samtools flagstat \
    "${PREFIX}_sorted.bam" \
    > "${PREFIX}.flagstat.txt"


echo "Done."
echo "BAM: ${OUTDIR}/${PREFIX}_sorted.bam"
echo "CSI: ${OUTDIR}/${PREFIX}_sorted.bam.csi"
