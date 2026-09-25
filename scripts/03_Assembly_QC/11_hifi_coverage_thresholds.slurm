#!/bin/bash
#SBATCH --job-name=hifi_cov_thresholds
#SBATCH --output=hifi_cov_thresholds_%j.out
#SBATCH --error=hifi_cov_thresholds_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

# Required software:
#   samtools
#   mosdepth

# -----------------------------
# Inputs
# -----------------------------

REF="/path/to/assembly_split.fa"
BAM="/path/to/hifi_cov_split_sorted.bam"
OUTDIR="/path/to/output"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

BED="${OUTDIR}/final_assembly_contigs.bed"
PREFIX="${OUTDIR}/hifi_thresholds"
SUMMARY="${OUTDIR}/hifi_thresholds.summary.tsv"

mkdir -p "${OUTDIR}"

# -----------------------------
# Create whole-contig regions
# -----------------------------

samtools faidx "${REF}"

awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' \
    "${REF}.fai" \
    > "${BED}"

# -----------------------------
# Calculate coverage thresholds
# -----------------------------

mosdepth \
    --threads "${THREADS}" \
    --no-per-base \
    --by "${BED}" \
    --thresholds 1,5,10,20,50 \
    "${PREFIX}" \
    "${BAM}"

# -----------------------------
# Summarize assembly-base coverage
# -----------------------------

zcat "${PREFIX}.thresholds.bed.gz" \
| awk '
BEGIN {
    OFS = "\t"
}

/^#/ {
    next
}

{
    total += $3 - $2
    x1  += $5
    x5  += $6
    x10 += $7
    x20 += $8
    x50 += $9
}

END {
    print "Threshold", "Bases", "Percent_assembly"
    printf ">=1x\t%d\t%.4f\n",  x1,  100*x1/total
    printf ">=5x\t%d\t%.4f\n",  x5,  100*x5/total
    printf ">=10x\t%d\t%.4f\n", x10, 100*x10/total
    printf ">=20x\t%d\t%.4f\n", x20, 100*x20/total
    printf ">=50x\t%d\t%.4f\n", x50, 100*x50/total
}
' > "${SUMMARY}"

echo "Coverage threshold summary: ${SUMMARY}"
