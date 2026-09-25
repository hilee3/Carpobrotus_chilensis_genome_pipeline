#!/bin/bash
#SBATCH --job-name=tiara_screening
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=tiara.%j.out
#SBATCH --error=tiara.%j.err

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Screen the chloroplast-filtered nuclear assembly using Tiara
# and remove residual non-nuclear contigs after manual validation.
#
# Tiara identified four mitochondrial contigs and one bacterial
# contig. These contigs were subsequently evaluated using
# sequencing depth and sequence-similarity evidence and removed
# from the nuclear assembly.

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
ASSEMBLY="/path/to/cp_cov_filter.filtered_nuclear.fa"
OUTDIR="tiara_output"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUTDIR}"

# -----------------------------
# Software requirements
# -----------------------------
# Tiara
# seqkit

# -----------------------------
# Check dependencies
# -----------------------------

for cmd in tiara seqkit; do
    command -v "${cmd}" >/dev/null || {
        echo "[ERROR] ${cmd} not found."
        exit 127
    }
done

[[ -s "${ASSEMBLY}" ]] || {
    echo "[ERROR] Assembly not found: ${ASSEMBLY}"
    exit 2
}

# -----------------------------
# STEP 1. Run Tiara
# -----------------------------

echo "[STEP 1] Running Tiara"

tiara \
    -i "${ASSEMBLY}" \
    -o "${OUTDIR}" \
    -t "${THREADS}" \
    --pr

# -----------------------------
# STEP 2. Remove validated non-nuclear contigs
# -----------------------------
#
# Tiara classifications were manually evaluated using read depth
# and sequence-similarity evidence. The following five contigs
# were removed from the final nuclear assembly:
#
#   ptg000045l_1  - mitochondrial
#   ptg000068c_1  - mitochondrial
#   ptg001036l_1  - mitochondrial
#   ptg001335l_1  - mitochondrial
#   ptg001958l_1  - bacterial

echo "[STEP 2] Removing validated non-nuclear contigs"

cat > tiara_remove_ids.txt << 'EOF'
ptg000045l_1
ptg000068c_1
ptg001036l_1
ptg001335l_1
ptg001958l_1
EOF

seqkit grep \
    -v \
    -f tiara_remove_ids.txt \
    "${ASSEMBLY}" \
    > cp_cov_filter.filtered_nuclear.tiara_filtered.fa

echo "[DONE] Tiara screening and filtering completed."
echo "[DONE] Removed contigs: $(wc -l < tiara_remove_ids.txt)"
echo "[DONE] Output: cp_cov_filter.filtered_nuclear.tiara_filtered.fa"