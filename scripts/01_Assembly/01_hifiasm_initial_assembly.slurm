#!/usr/bin/env bash
#SBATCH --job-name=hifiasm_chilensis
#SBATCH --time=10-01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=hifiasm_chilensis.out
#SBATCH --error=hifiasm_chilensis.err

set -euo pipefail

# Generate a de novo genome assembly of Carpobrotus chilensis
# from PacBio HiFi reads using hifiasm.

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
READS="/path/to/pacbio_hifi_reads.fastq.gz"
PREFIX="chilensis_assembly_l0_031025.asm"

THREADS="${SLURM_CPUS_PER_TASK:-32}"

# -----------------------------
# Software requirements
# -----------------------------
# hifiasm must be available in PATH.

hifiasm \
    -o "${PREFIX}" \
    -t "${THREADS}" \
    "${READS}"