#!/usr/bin/env bash
#SBATCH --job-name=hifi_qc
#SBATCH --output=hifi_qc_%j.out
#SBATCH --error=hifi_qc_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# Required software:
#   samtools
#   mosdepth

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/output"
BAM="${WORKDIR}/hifi_cov_split_sorted.bam"

FLAGSTAT="${WORKDIR}/hifi_cov_split.flagstat.txt"
STATS="${WORKDIR}/hifi_cov_split.stats.txt"

MOSPREFIX="${WORKDIR}/hifi_cov_1kb"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

cd "${WORKDIR}"

# -----------------------------
# Check input
# -----------------------------

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}" >&2
    exit 1
fi

# -----------------------------
# Check BAM
# -----------------------------

samtools quickcheck "${BAM}"

# -----------------------------
# 1. Read mapping statistics
# -----------------------------

samtools flagstat \
    -@ "${THREADS}" \
    "${BAM}" \
    > "${FLAGSTAT}"

# -----------------------------
# 2. Detailed alignment statistics
# -----------------------------

samtools stats \
    "${BAM}" \
    > "${STATS}"

# -----------------------------
# 3. Depth distribution in 1-kb windows
# -----------------------------

mosdepth \
    --threads "${THREADS}" \
    --no-per-base \
    --by 1000 \
    "${MOSPREFIX}" \
    "${BAM}"

echo "Done."
echo "Flagstat: ${FLAGSTAT}"
echo "Samtools stats: ${STATS}"
echo "Mosdepth summary: ${MOSPREFIX}.mosdepth.summary.txt"
echo "1-kb depth: ${MOSPREFIX}.regions.bed.gz"
echo "Depth distribution: ${MOSPREFIX}.mosdepth.global.dist.txt"
