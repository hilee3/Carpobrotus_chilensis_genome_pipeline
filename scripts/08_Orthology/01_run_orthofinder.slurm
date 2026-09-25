#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --output=orthofinder_%j.out
#SBATCH --error=orthofinder_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -euo pipefail

# Identify orthogroups between C. chilensis and C. edulis
# using the longest predicted protein isoform per gene.

CHILENSIS_PROTEIN="/path/to/braker.g103t2_removed.longest_isoform.gffread.aa"
EDULIS_PROTEIN="/path/to/braker.longest_isoform.gffread.aa"

WORK_DIR="/path/to/genome_comparison"
PROTEOME_DIR="${WORK_DIR}/proteomes"

mkdir -p "${PROTEOME_DIR}"

cp "${CHILENSIS_PROTEIN}" "${PROTEOME_DIR}/chilensis.fa"
cp "${EDULIS_PROTEIN}" "${PROTEOME_DIR}/edulis.fa"

cd "${WORK_DIR}"

orthofinder \
    -f "${PROTEOME_DIR}" \
    -t 8 \
    -a 2 \
    -S diamond \
    -o orthofinder_results

echo "OrthoFinder analysis completed"
echo "Results in: ${WORK_DIR}/orthofinder_results"