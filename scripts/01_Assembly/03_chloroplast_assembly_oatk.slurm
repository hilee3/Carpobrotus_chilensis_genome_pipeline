#!/usr/bin/env bash
#SBATCH --job-name=oatk_plastid
#SBATCH --output=oatk_plastid_%j.out
#SBATCH --error=oatk_plastid_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

# Assemble the C. chilensis chloroplast genome from PacBio HiFi reads
# using Oatk.

# ===== Paths =====
READS="/path/to/pacbio_hifi_reads.fastq.gz"              # HiFi read file
PLASTID_DB="/path/to/oatkdb/embryophyta_pltd.fam"        # Chloroplast HMM DB
OUTPUT_PREFIX="carp_plastid"                              # Output file prefix
NHMM_PATH=$(which nhmmscan)

# ===== Parameters =====
KMER=1001
CUTOFF=400
THREADS=16

# ===== Run Oatk =====
echo "Starting Oatk chloroplast assembly..."

oatk \
  -k "${KMER}" \
  -c "${CUTOFF}" \
  -t "${THREADS}" \
  --nhmmscan "${NHMM_PATH}" \
  -p "${PLASTID_DB}" \
  -o "${OUTPUT_PREFIX}" \
  "${READS}"

echo "Assembly finished. Main outputs:"
echo " - ${OUTPUT_PREFIX}.utg.final.gfa"
echo " - ${OUTPUT_PREFIX}.pltd.ctg.fasta (final chloroplast contig)"
echo " - ${OUTPUT_PREFIX}.pltd.bed / .gfa (annotations & structure)"
