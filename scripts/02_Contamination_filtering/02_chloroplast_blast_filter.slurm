#!/usr/bin/env bash
#SBATCH --job-name=chloroplast_blast_filter
#SBATCH --output=chloroplast_blast_filter_%j.out
#SBATCH --error=chloroplast_blast_filter_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Remove residual chloroplast-like contigs using BLASTN.
#
# Contigs are removed if they satisfy:
#   - percent identity >= 99
#   - query coverage >= 99
#
# Input:
#   - organelle-filtered nuclear assembly (nuclear_cleaned_v2.fa)
#   - assembled chloroplast genome
#
# Output:
#   - BLASTN table
#   - list of removed contigs
#   - filtered nuclear assembly

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
CP_FASTA="/path/to/assembled_chloroplast_genome.fa"
NUC_FASTA="/path/to/nuclear_cleaned_v2.fa"

PREFIX="cp_cov_filter"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# -----------------------------
# Software requirements
# -----------------------------
# BLAST+
# seqkit

# -----------------------------
# Check dependencies
# -----------------------------
for cmd in makeblastdb blastn seqkit awk sort; do
    command -v "${cmd}" >/dev/null || {
        echo "[ERROR] ${cmd} not found."
        exit 127
    }
done

[[ -s "${CP_FASTA}" ]] || { echo "[ERROR] CP_FASTA not found: ${CP_FASTA}"; exit 2; }
[[ -s "${NUC_FASTA}" ]] || { echo "[ERROR] NUC_FASTA not found: ${NUC_FASTA}"; exit 2; }

echo "[INFO] Chloroplast genome: ${CP_FASTA}"
echo "[INFO] Nuclear assembly: ${NUC_FASTA}"

echo "[STEP 1] Building chloroplast BLAST database"
makeblastdb -in "${CP_FASTA}" -dbtype nucl -out "${PREFIX}.cpdb"

echo "[STEP 2] Running BLASTN"
blastn \
    -query "${NUC_FASTA}" \
    -db "${PREFIX}.cpdb" \
    -out "${PREFIX}.blastn.tsv" \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovs" \
    -evalue 1e-20 \
    -num_threads "${THREADS}"

echo "[STEP 3] Identifying high-confidence chloroplast-like contigs"
awk 'BEGIN{FS=OFS="\t"}
$3 >= 99 && $12 >= 99 {print}' \
    "${PREFIX}.blastn.tsv" \
    > "${PREFIX}.blastn.filtered.tsv"

awk '{print $1}' "${PREFIX}.blastn.filtered.tsv" \
    | sort -u \
    > "${PREFIX}.contigs_to_remove.txt"

echo "[INFO] Contigs to remove: $(wc -l < "${PREFIX}.contigs_to_remove.txt")"

echo "[STEP 4] Removing chloroplast-like contigs"
seqkit grep \
    -v \
    -f "${PREFIX}.contigs_to_remove.txt" \
    "${NUC_FASTA}" \
    > "${PREFIX}.filtered_nuclear.fa"

echo "[DONE] Final output assembly: ${PREFIX}.filtered_nuclear.fa"