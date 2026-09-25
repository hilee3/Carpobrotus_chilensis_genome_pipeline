#!/usr/bin/env bash
#SBATCH --job-name=braker_edulis_etp
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --output=braker_edulis_%j.out
#SBATCH --error=braker_edulis_%j.err

set -euo pipefail

# Required software:
#   BRAKER3
#   AUGUSTUS
#   GeneMark-ETP
#   ProtHint
#   DIAMOND
#   StringTie
#   samtools
#   bamtools
#   bedtools
#   gffread

# -----------------------------
# Inputs
# -----------------------------

GENOME="/path/to/edulis_nuclear_only.softmasked.fa"

BAM1="/path/to/ERR15315142_1.sorted.bam"
BAM2="/path/to/ERR15315144_1.sorted.bam"
BAM="${BAM1},${BAM2}"

# Viridiplantae protein evidence.
PROT_SEQ="/path/to/Viridiplantae.fa"

FINAL_OUT="/path/to/output/edulis_braker3"

SPECIES="carpobrotus_edulis_etp"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# -----------------------------
# Software-specific paths
# -----------------------------

GENEMARK_PATH="/path/to/gmetp_linux_64/bin"
PROTHINT_PATH="/path/to/ProtHint/bin"

export GENEMARK_PATH
export PROTHINT_PATH
export PATH="${GENEMARK_PATH}:${PROTHINT_PATH}:${PATH}"

# -----------------------------
# Scratch workspace
# -----------------------------

SCRATCH_BASE="/path/to/scratch"
JOBTMP="${SCRATCH_BASE}/braker_${SLURM_JOB_ID:-$$}"
WORK="${JOBTMP}/work"

mkdir -p "${WORK}" "${FINAL_OUT}"

export TMPDIR="${JOBTMP}"
export TEMP="${JOBTMP}"
export TMP="${JOBTMP}"

# -----------------------------
# Run BRAKER3
# -----------------------------

braker.pl \
    --genome="${GENOME}" \
    --prot_seq="${PROT_SEQ}" \
    --bam="${BAM}" \
    --species="${SPECIES}" \
    --softmasking \
    --threads="${THREADS}" \
    --workingdir="${WORK}" \
    --GENEMARK_PATH="${GENEMARK_PATH}" \
    --gff3

# -----------------------------
# Copy primary outputs
# -----------------------------

for f in \
    braker.gtf \
    braker.gff3 \
    braker.aa \
    braker.codingseq \
    hintsfile.gff \
    what-to-cite.txt \
    braker.log
do
    [[ -e "${WORK}/${f}" ]] && cp -a "${WORK}/${f}" "${FINAL_OUT}/"
done

echo "[INFO] Output: ${FINAL_OUT}"
