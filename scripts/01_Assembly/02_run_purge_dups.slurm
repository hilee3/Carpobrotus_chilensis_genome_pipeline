#!/usr/bin/env bash
#SBATCH --job-name=purge_dups
#SBATCH --output=purge_dups.%j.out
#SBATCH --error=purge_dups.%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

set -euo pipefail

# Remove residual haplotypic redundancy from the hifiasm primary
# assembly using assembly self-alignment.
#
# Coverage-based sequence classification was not used.

ASM="/path/to/hifiasm_primary_assembly.fa"
PREFIX="pb.083025"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# Split the primary assembly.
split_fa "${ASM}" > "${PREFIX}.split.fa"

# Self-align the primary assembly.
minimap2 \
    -x asm5 \
    -DP \
    -I 50G \
    -t "${THREADS}" \
    "${PREFIX}.split.fa" \
    "${PREFIX}.split.fa" \
    > "${PREFIX}.split.self.paf"

# Identify duplicated regions.
purge_dups \
    -2 \
    "${PREFIX}.split.self.paf" \
    > dups.bed \
    2> purge_dups.log

# Extract the purged primary assembly and removed haplotig sequences.
get_seqs \
    -e \
    dups.bed \
    "${ASM}"

mv purged.fa purged_primary.fa
mv hap.fa hap_from_primary.fa  # Public release: Cchil_purge_dups_removed_haplotigs.fa