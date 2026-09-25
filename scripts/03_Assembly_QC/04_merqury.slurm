#!/usr/bin/env bash
#SBATCH --job-name=merqury
#SBATCH --output=merqury.%j.out
#SBATCH --error=merqury.%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

set -euo pipefail

# Assess consensus quality and k-mer completeness of the
# final genome assembly using Merqury.

# -----------------------------
# Inputs
# -----------------------------
MERYL_DB="/path/to/pacbio.meryl"
ASSEMBLY="/path/to/assembly_split.fa"

OUT_PREFIX="merqury"
WORKDIR="merqury_results"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

export OMP_NUM_THREADS="${THREADS}"

# -----------------------------
# Software requirements
# -----------------------------
# Merqury
# meryl

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

meryl statistics "${MERYL_DB}" \
    | tee "${OUT_PREFIX}.db.stats.txt"

meryl histogram "${MERYL_DB}" \
    > "${OUT_PREFIX}.db.hist"

merqury.sh \
    "${MERYL_DB}" \
    "${ASSEMBLY}" \
    "${OUT_PREFIX}"

echo "[DONE] Merqury analysis completed."