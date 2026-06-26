#!/usr/bin/env bash
#SBATCH --job-name=run_purge_dups
#SBATCH --output=run_purge_dups.%j.out
#SBATCH --error=run_purge_dups.%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Remove duplicated sequence from the hifiasm primary assembly
# using purge_dups.
#
# Workflow:
#   1. Split the primary assembly.
#   2. Self-align the primary assembly.
#   3. Run purge_dups on the primary assembly.
#   4. Extract the purged primary assembly and haplotig sequences.
#   5. Merge the extracted haplotigs with the hifiasm hap1 and
#      hap2 assemblies.
#   6. Perform a second purge_dups run on the merged haplotigs.
#
# Coverage cutoffs were initially generated using calcuts and
# subsequently adjusted after manual inspection of the PB.stat
# coverage histogram. The final cutoff values used in this study
# were:
#
#     30 80 250


# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
ASM="/path/to/hifiasm_primary_assembly.fa"
HAP1="/path/to/hifiasm_hap1_assembly.fa"
HAP2="/path/to/hifiasm_hap2_assembly.fa"

# Coverage statistics and manually adjusted cutoff file.
COV="PB.base.cov"
CUT="cutoffs_custom"

PREFIX="pb.083025"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# -----------------------------
# Software requirements
# -----------------------------
# purge_dups utilities: purge_dups, split_fa, get_seqs
# minimap2

# -----------------------------
# Check dependencies
# -----------------------------

command -v purge_dups >/dev/null || { echo "[ERROR] purge_dups not found in PATH"; exit 127; }
command -v minimap2   >/dev/null || { echo "[ERROR] minimap2 not found in PATH"; exit 127; }
command -v split_fa   >/dev/null || { echo "[ERROR] split_fa not found in PATH"; exit 127; }
command -v get_seqs   >/dev/null || { echo "[ERROR] get_seqs not found in PATH"; exit 127; }

# -----------------------------
# Check inputs
# -----------------------------

[[ -s "${ASM}"  ]] || { echo "[ERROR] Assembly not found: ${ASM}"; exit 2; }
[[ -s "${HAP1}" ]] || { echo "[ERROR] hap1 assembly not found: ${HAP1}"; exit 2; }
[[ -s "${HAP2}" ]] || { echo "[ERROR] hap2 assembly not found: ${HAP2}"; exit 2; }
[[ -s "${COV}"  ]] || { echo "[ERROR] Coverage file not found: ${COV}"; exit 2; }
[[ -s "${CUT}"  ]] || { echo "[ERROR] Cutoff file not found: ${CUT}"; exit 2; }

echo "[INFO] Primary assembly : ${ASM}"
echo "[INFO] Coverage file    : ${COV}"
echo "[INFO] Cutoff file      : ${CUT}"
echo "[INFO] Threads          : ${THREADS}"

# -----------------------------
# STEP 1-4: Purge primary assembly
# -----------------------------

echo "[STEP 1] Splitting primary assembly"
split_fa "${ASM}" > "${PREFIX}.split.fa"

echo "[STEP 2] Self-aligning primary assembly"
minimap2 \
    -x asm5 \
    -DP \
    -I 50G \
    -t "${THREADS}" \
    "${PREFIX}.split.fa" \
    "${PREFIX}.split.fa" \
    > "${PREFIX}.split.self.paf"

echo "[STEP 3] Running purge_dups on the primary assembly"
purge_dups \
    -2 \
    -T "${CUT}" \
    -c "${COV}" \
    "${PREFIX}.split.self.paf" \
    > dups.bed \
    2> purge_dups.log

echo "[STEP 4] Extracting purged primary assembly and haplotigs"
get_seqs -e dups.bed "${ASM}"

mv purged.fa purged_primary.fa
mv hap.fa hap_from_primary.fa

# -----------------------------
# STEP 5-9: Refine haplotigs
# -----------------------------

echo "[STEP 5] Merging haplotigs with hifiasm hap1 and hap2 assemblies"
cat \
    hap_from_primary.fa \
    "${HAP1}" \
    "${HAP2}" \
    > merged_hap.fa

echo "[STEP 6] Splitting merged haplotigs"
split_fa merged_hap.fa > merged_hap.split.fa

echo "[STEP 7] Self-aligning merged haplotigs"
minimap2 \
    -x asm5 \
    -DP \
    -I 50G \
    -t "${THREADS}" \
    merged_hap.split.fa \
    merged_hap.split.fa \
    > merged_hap.self.paf

echo "[STEP 8] Running purge_dups on merged haplotigs"
purge_dups \
    -T "${CUT}" \
    -c "${COV}" \
    merged_hap.self.paf \
    > hap_dups.bed \
    2> hap_purge.log

echo "[STEP 9] Extracting refined haplotigs"
get_seqs -e hap_dups.bed merged_hap.fa

mv purged.fa refined_haplotigs.fa
mv hap.fa filtered_out_from_merged.fa

echo "[DONE] Purged primary assembly        : purged_primary.fa"
echo "[DONE] Refined haplotigs              : refined_haplotigs.fa"
echo "[DONE] Filtered haplotigs             : filtered_out_from_merged.fa"