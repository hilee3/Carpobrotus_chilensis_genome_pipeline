#!/bin/bash
#SBATCH --job-name=TE_gene_overlap
#SBATCH --output=TE_gene_overlap_%j.out
#SBATCH --error=TE_gene_overlap_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

# Required software:
#   bedtools

# -----------------------------
# Inputs
# -----------------------------

RM_OUT="/path/to/genome.fa.out"   # RepeatMasker output from the final repeat annotation.
GENE_GFF="/path/to/braker.g103t2_removed.longest_isoform.agat.gff3"
OUTDIR="/path/to/output/TE_gene_overlap"

mkdir -p "${OUTDIR}"

TE_BED="${OUTDIR}/repeatmasker_TE.bed"
GENE_BED="${OUTDIR}/genes.bed"
CDS_BED="${OUTDIR}/CDS.bed"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools not found in PATH"
    exit 1
}

# ------------------------------------------------------------
# 1. Inspect RepeatMasker repeat classes
# ------------------------------------------------------------

awk '
NR > 3 && NF >= 11 {
    print $11
}' "${RM_OUT}" \
    | sort \
    | uniq -c \
    | sort -nr \
    > "${OUTDIR}/repeatmasker_classes.txt"

# ------------------------------------------------------------
# 2. Extract TE intervals from RepeatMasker .out
#
# Include:
#   LTR
#   LINE
#   SINE
#   DNA transposons
#   RC/Helitron
#   Unknown/unclassified interspersed repeats
#
# Exclude:
#   Simple_repeat
#   Low_complexity
#
# BED uses 0-based start coordinates.
# ------------------------------------------------------------

awk 'BEGIN{OFS="\t"}
NR > 3 && NF >= 11 {
    class=$11

    if (class != "Simple_repeat" && class != "Low_complexity") {
        print $5, $6-1, $7, $10, class
    }
}' "${RM_OUT}" \
    | sort -k1,1 -k2,2n \
    > "${TE_BED}"

# ------------------------------------------------------------
# 3. Extract genes from final longest-isoform GFF3
# ------------------------------------------------------------

awk -F'\t' '
BEGIN {OFS="\t"}
$0 !~ /^#/ && $3=="gene" {
    id="."
    n=split($9,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^ID=/){
            sub(/^ID=/,"",a[i])
            id=a[i]
        }
    }
    print $1,$4-1,$5,id
}' "${GENE_GFF}" \
    | sort -k1,1 -k2,2n \
    > "${GENE_BED}"

# ------------------------------------------------------------
# 4. Extract CDS intervals
#    Keep Parent transcript ID for gene-level grouping.
# ------------------------------------------------------------

awk -F'\t' '
BEGIN {OFS="\t"}
$0 !~ /^#/ && $3=="CDS" {
    parent="."
    n=split($9,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^Parent=/){
            sub(/^Parent=/,"",a[i])
            parent=a[i]
        }
    }
    print $1,$4-1,$5,parent
}' "${GENE_GFF}" \
    | sort -k1,1 -k2,2n \
    > "${CDS_BED}"

# ------------------------------------------------------------
# 5. Genes with any overlap with a TE
# ------------------------------------------------------------

bedtools intersect \
    -a "${GENE_BED}" \
    -b "${TE_BED}" \
    -u \
    > "${OUTDIR}/genes_overlapping_TE.bed"

TOTAL_GENES=$(wc -l < "${GENE_BED}")
GENES_TE=$(wc -l < "${OUTDIR}/genes_overlapping_TE.bed")

# ------------------------------------------------------------
# 6. CDS intervals overlapping TEs
# ------------------------------------------------------------

bedtools intersect \
    -a "${CDS_BED}" \
    -b "${TE_BED}" \
    -u \
    > "${OUTDIR}/CDS_overlapping_TE.bed"

cut -f4 "${OUTDIR}/CDS_overlapping_TE.bed" \
    | sort -u \
    > "${OUTDIR}/transcripts_with_CDS_TE_overlap.txt"

CDS_TE_GENES=$(
    wc -l < "${OUTDIR}/transcripts_with_CDS_TE_overlap.txt"
)

# ------------------------------------------------------------
# 7. Fraction of CDS sequence overlapping TEs per transcript
#
# Merge TE intervals first to avoid double-counting
# overlapping TE annotations.
# ------------------------------------------------------------

bedtools merge \
    -i "${TE_BED}" \
    > "${OUTDIR}/repeatmasker_TE.merged.bed"

# Total CDS length per transcript
awk '
{
    len=$3-$2
    total[$4]+=len
}
END {
    for(id in total)
        print id,total[id]
}' OFS="\t" "${CDS_BED}" \
    > "${OUTDIR}/CDS_total_length.tsv"

# Number of CDS bp overlapping merged TE intervals
bedtools intersect \
    -a "${CDS_BED}" \
    -b "${OUTDIR}/repeatmasker_TE.merged.bed" \
    -wo \
    > "${OUTDIR}/CDS_TE_intersections.tsv"

awk '
{
    id=$4
    overlap[id]+=$NF
}
END {
    for(id in overlap)
        print id,overlap[id]
}' OFS="\t" "${OUTDIR}/CDS_TE_intersections.tsv" \
    > "${OUTDIR}/CDS_TE_overlap_bp.tsv"

# Combine total CDS length and TE-overlap length
awk '
BEGIN {OFS="\t"}

FNR==NR {
    total[$1]=$2
    next
}

{
    overlap[$1]=$2
}

END {
    print "ID","CDS_length_bp","TE_overlap_bp","TE_overlap_percent"

    for(id in total) {
        ov=(id in overlap ? overlap[id] : 0)
        pct=100*ov/total[id]
        print id,total[id],ov,pct
    }
}' "${OUTDIR}/CDS_total_length.tsv" \
   "${OUTDIR}/CDS_TE_overlap_bp.tsv" \
   > "${OUTDIR}/CDS_TE_overlap_summary.tsv"

# ------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------

echo "========================================"
echo "SUMMARY"
echo "========================================"

printf "Total predicted genes: %d\n" "${TOTAL_GENES}"

awk -v x="${GENES_TE}" -v n="${TOTAL_GENES}" \
    'BEGIN{
        printf "Genes overlapping any TE: %d (%.2f%%)\n", x,100*x/n
    }'

awk -v x="${CDS_TE_GENES}" -v n="${TOTAL_GENES}" \
    'BEGIN{
        printf "Genes with CDS overlapping TE: %d (%.2f%%)\n", x,100*x/n
    }'

echo
echo "CDS TE-overlap categories:"

awk -F'\t' '
NR>1 {
    n++

    if($4 == 0)
        zero++
    else if($4 <= 10)
        low++
    else if($4 <= 50)
        medium++
    else
        high++
}
END {
    printf "0%% overlap:       %d (%.2f%%)\n", zero,100*zero/n
    printf ">0-10%% overlap:   %d (%.2f%%)\n", low,100*low/n
    printf ">10-50%% overlap:  %d (%.2f%%)\n", medium,100*medium/n
    printf ">50%% overlap:     %d (%.2f%%)\n", high,100*high/n
}' "${OUTDIR}/CDS_TE_overlap_summary.tsv"
