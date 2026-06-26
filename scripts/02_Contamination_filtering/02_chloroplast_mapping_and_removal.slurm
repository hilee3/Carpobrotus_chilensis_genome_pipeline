#!/usr/bin/env bash
#SBATCH --job-name=cp_mapping_removal
#SBATCH --output=cp_mapping_removal.%j.out
#SBATCH --error=cp_mapping_removal.%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Screen the organelle-filtered nuclear assembly for residual
# chloroplast-derived contigs and remove high-confidence
# chloroplast-like contigs.
#
# This script aligns the nuclear assembly against the assembled
# chloroplast genome, summarizes alignment coverage for each
# contig, and removes contigs with best alignment coverage >= 0.80.
#
# Input:
#   - organelle-filtered nuclear assembly (nuclear_cleaned_v2.fa)
#   - assembled chloroplast genome
#
# Output:
#   - alignment summary tables
#   - remove_cp_ids.txt
#   - nuclear_cleaned_v3.fa

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
QUERY="/path/to/nuclear_cleaned_v2.fa"
CP_REF="/path/to/chloroplast_genome.fa"

OUTDIR="cp_mapping_check"

PREFIX="nuclear_vs_cp"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# -----------------------------
# Software requirements
# -----------------------------
# minimap2
# samtools
# seqkit

# -----------------------------
# Check dependencies
# -----------------------------

for cmd in minimap2 samtools awk sort seqkit; do
    command -v "${cmd}" >/dev/null || {
        echo "[ERROR] ${cmd} not found."
        exit 127
    }
done

echo "[INFO] Query assembly: ${QUERY}"
echo "[INFO] Chloroplast genome: ${CP_REF}"

# -----------------------------
# STEP 1. Index chloroplast genome
# -----------------------------

echo "[STEP 1] Indexing chloroplast genome"

samtools faidx "${CP_REF}"

# -----------------------------
# STEP 2. Align nuclear contigs
# -----------------------------

echo "[STEP 2] Aligning nuclear contigs"

minimap2 \
    -x asm5 \
    -c \
    -t "${THREADS}" \
    "${CP_REF}" \
    "${QUERY}" \
    > "${PREFIX}.paf"

minimap2 \
    -ax asm5 \
    -t "${THREADS}" \
    "${CP_REF}" \
    "${QUERY}" \
| samtools sort \
    -@ "${THREADS}" \
    -o "${PREFIX}.bam" -

samtools index "${PREFIX}.bam"

# -----------------------------
# STEP 3. Extract contig lengths
# -----------------------------

echo "[STEP 3] Extracting contig lengths"

samtools faidx "${QUERY}"
cut -f1,2 "${QUERY}.fai" > query_lengths.tsv

# -----------------------------
# STEP 4. Summarize alignments
# -----------------------------

echo "[STEP 4] Summarizing alignments"

awk 'BEGIN{
    OFS="\t";
    print "contig","qlen","sum_aln_len","best_aln_len","best_mapq","sum_cov_ratio","best_cov_ratio","target_names"
}
{
    q=$1; qlen=$2; t=$6; aln=$11; mapq=$12;

    qlen_map[q]=qlen;
    sum_aln[q]+=aln;

    if(aln>best_aln[q]) best_aln[q]=aln;
    if(mapq>best_mapq[q]) best_mapq[q]=mapq;

    if(!(seen[q,t]++)){
        if(targets[q]=="") targets[q]=t;
        else targets[q]=targets[q]","t;
    }
}
END{
    for(q in qlen_map){
        sum_cov=(qlen_map[q]>0?sum_aln[q]/qlen_map[q]:0);
        best_cov=(qlen_map[q]>0?best_aln[q]/qlen_map[q]:0);

        print q,
              qlen_map[q],
              sum_aln[q]+0,
              best_aln[q]+0,
              best_mapq[q]+0,
              sum_cov,
              best_cov,
              targets[q];
    }
}' "${PREFIX}.paf" \
| sort -k6,6nr \
> "${PREFIX}.alignment_summary.tsv"

# -----------------------------
# STEP 5. Add contigs without hits
# -----------------------------

echo "[STEP 5] Adding contigs without chloroplast hits"

awk 'NR==FNR {hit[$1]=1; next}
     !($1 in hit){
         print $1"\t"$2"\t0\t0\t0\t0\t0\tNA
     }' \
    "${PREFIX}.alignment_summary.tsv" \
    query_lengths.tsv \
    > "${PREFIX}.nohit.tsv"

{
    head -n1 "${PREFIX}.alignment_summary.tsv"
    tail -n+2 "${PREFIX}.alignment_summary.tsv"
    cat "${PREFIX}.nohit.tsv"
} > "${PREFIX}.all_contigs.tsv"

{
    head -n1 "${PREFIX}.all_contigs.tsv"
    tail -n+2 "${PREFIX}.all_contigs.tsv" | sort -k6,6nr
} > "${PREFIX}.all_contigs.sorted.tsv"

# -----------------------------
# STEP 6. Candidate tables
# -----------------------------

echo "[STEP 6] Generating candidate tables"

awk 'BEGIN{FS=OFS="\t"}
NR==1{print;next}
$6>=0.80 || $7>=0.80{print}' \
"${PREFIX}.all_contigs.sorted.tsv" \
> "${PREFIX}.strong_cp_hits.tsv"

awk 'BEGIN{FS=OFS="\t"}
NR==1{print;next}
($6>=0.20 && $6<0.80) || ($7>=0.20 && $7<0.80){print}' \
"${PREFIX}.all_contigs.sorted.tsv" \
> "${PREFIX}.partial_cp_hits.tsv"

awk 'BEGIN{FS=OFS="\t"}
NR==1{print;next}
$6<0.20 && $7<0.20{print}' \
"${PREFIX}.all_contigs.sorted.tsv" \
> "${PREFIX}.weak_or_no_cp_hits.tsv"

# -----------------------------
# STEP 7. Summary
# -----------------------------

echo "[STEP 7] Summary"

TOTAL=$(awk 'END{print NR-1}' "${PREFIX}.all_contigs.sorted.tsv")
STRONG=$(awk 'END{print NR-1}' "${PREFIX}.strong_cp_hits.tsv")
PARTIAL=$(awk 'END{print NR-1}' "${PREFIX}.partial_cp_hits.tsv")
WEAK=$(awk 'END{print NR-1}' "${PREFIX}.weak_or_no_cp_hits.tsv")

cat <<EOF

=========================================
Chloroplast mapping summary
=========================================

Total contigs:            ${TOTAL}
Strong CP hits:           ${STRONG}
Partial CP hits:          ${PARTIAL}
Weak or no CP hits:       ${WEAK}

Output tables:

${PREFIX}.all_contigs.sorted.tsv
${PREFIX}.strong_cp_hits.tsv
${PREFIX}.partial_cp_hits.tsv
${PREFIX}.weak_or_no_cp_hits.tsv

=========================================

EOF

# -----------------------------
# STEP 8. Remove high-confidence chloroplast-like contigs
# -----------------------------
echo "[STEP 8] Removing high-confidence chloroplast-like contigs"

awk 'BEGIN{FS=OFS="\t"}
NR>1 && $7 >= 0.80 {print $1}' \
    "${PREFIX}.all_contigs.sorted.tsv" \
    > remove_cp_ids.txt

seqkit grep \
    -v \
    -f remove_cp_ids.txt \
    "${QUERY}" \
    > nuclear_cleaned_v3.fa

echo "[DONE] Removed chloroplast-like contigs: $(wc -l < remove_cp_ids.txt)"
echo "[DONE] Final output assembly: nuclear_cleaned_v3.fa"
