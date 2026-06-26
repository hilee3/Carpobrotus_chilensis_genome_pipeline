#!/usr/bin/env bash
#SBATCH --job-name=organelle_filtering
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=organelle_filtering_%j.out
#SBATCH --error=organelle_filtering_%j.err

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Identify and remove putative plastid- and mitochondrial-
# derived contigs from the purge_dups-filtered primary assembly.
#
# Criteria:
#   1. Alignment coverage to plastid or mitochondrial reference >= 0.50
#   2. Per-contig sequencing depth >= 3x the nuclear median depth
#   3. Sequencing-depth ratio relative to organelle references:
#        plastid >= 0.10
#        mitochondrion >= 0.05
#
# Organelle references:
#   Plastid and mitochondrial sequences from the Carpobrotus edulis
#   reference assembly, NCBI assembly GCA_965788445.1.

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
ASM="/path/to/purge_dups_filtered_primary_assembly.fa"
READS="/path/to/pacbio_hifi_reads.fastq.gz"

CHLORO="/path/to/C_edulis_plastid.fa"
MITO="/path/to/C_edulis_mitochondria.fa"

OUTDIR="organelle_filtering"

READS_TYPE="hifi"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# -----------------------------
# Filtering thresholds
# -----------------------------

ALCR_THRESHOLD="0.50"
R_NUC_MIN="3.0"

SDR_PLASTID="0.10"
SDR_MITO="0.05"

REMOVAL_MODE="remove"

# -----------------------------
# Software requirements
# -----------------------------
# minimap2
# samtools
# BLAST+
# seqtk
# python3

# -----------------------------
# Check dependencies
# -----------------------------

function need() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "[FATAL] $1 not found on PATH"
        exit 1
    }
}

need minimap2
need samtools
need blastn
need makeblastdb
need seqtk
need python3

# -----------------------------
# Prepare output directory
# -----------------------------

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "[INFO] Assembly: ${ASM}"
echo "[INFO] Reads: ${READS}"
echo "[INFO] Plastid reference: ${CHLORO}"
echo "[INFO] Mitochondrial reference: ${MITO}"
echo "[INFO] Threads: ${THREADS}"

[[ -s "${ASM}" ]] || { echo "[FATAL] ASM not found: ${ASM}"; exit 1; }
[[ -s "${READS}" ]] || { echo "[FATAL] READS not found: ${READS}"; exit 1; }
[[ -s "${CHLORO}" ]] || { echo "[FATAL] CHLORO not found: ${CHLORO}"; exit 1; }
[[ -s "${MITO}" ]] || { echo "[FATAL] MITO not found: ${MITO}"; exit 1; }

if [[ "${READS_TYPE}" == "hifi" ]]; then
    MAP_PRESET="-ax map-hifi"
elif [[ "${READS_TYPE}" == "illumina" ]]; then
    MAP_PRESET="-ax sr"
else
    echo "[FATAL] READS_TYPE must be hifi or illumina"
    exit 1
fi

# -----------------------------
# STEP 1: Depth on nuclear assembly
# -----------------------------

echo "[STEP 1] Calculating per-contig depth for the nuclear assembly"

minimap2 ${MAP_PRESET} -t "${THREADS}" "${ASM}" "${READS}" \
    | samtools sort -@ "${THREADS}" -o nuc.bam

samtools index nuc.bam

samtools depth -a nuc.bam \
    | awk '{sum[$1]+=$3; n[$1]++} END{for(c in sum) print c"\t"sum[c]/n[c]}' \
    > depth_nuc.tsv

python3 - << 'PY'
import statistics
import sys

vals = []
with open("depth_nuc.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            try:
                vals.append(float(parts[1]))
            except ValueError:
                pass

if not vals:
    print("[FATAL] No per-contig depths in depth_nuc.tsv", file=sys.stderr)
    sys.exit(1)

dnuc = statistics.median(vals)
open("D_nuc.txt", "w").write(f"{dnuc}\n")
print(f"[INFO] D_nuc = {dnuc:.3f}x")
PY

D_NUC=$(cat D_nuc.txt)

# -----------------------------
# STEP 2: Depth on organelle references
# -----------------------------

echo "[STEP 2] Calculating plastid and mitochondrial read depth"

minimap2 ${MAP_PRESET} -t "${THREADS}" "${CHLORO}" "${READS}" \
    | samtools sort -@ "${THREADS}" -o plastid.bam

samtools index plastid.bam

samtools depth -a plastid.bam \
    | awk '{sum+=$3; n++} END{print sum/n}' \
    > D_plastid.txt

minimap2 ${MAP_PRESET} -t "${THREADS}" "${MITO}" "${READS}" \
    | samtools sort -@ "${THREADS}" -o mito.bam

samtools index mito.bam

samtools depth -a mito.bam \
    | awk '{sum+=$3; n++} END{print sum/n}' \
    > D_mito.txt

D_PLASTID=$(cat D_plastid.txt)
D_MITO=$(cat D_mito.txt)

awk -v p="${D_PLASTID}" -v m="${D_MITO}" \
    'BEGIN{printf("D_plastid\t%.3f\nD_mito\t%.3f\n", p, m)}' \
    > depth_org_summary.txt

echo "[INFO] D_plastid = ${D_PLASTID}x"
echo "[INFO] D_mito = ${D_MITO}x"

# -----------------------------
# STEP 3: BLAST against organelle references
# -----------------------------

echo "[STEP 3] Running BLASTN against plastid and mitochondrial references"

makeblastdb -in "${CHLORO}" -dbtype nucl -out plastid_db

blastn -query "${ASM}" -db plastid_db \
    -out asm_vs_plastid.tsv \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend evalue bitscore' \
    -evalue 1e-20 \
    -num_threads "${THREADS}"

makeblastdb -in "${MITO}" -dbtype nucl -out mito_db

blastn -query "${ASM}" -db mito_db \
    -out asm_vs_mito.tsv \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend evalue bitscore' \
    -evalue 1e-20 \
    -num_threads "${THREADS}"

awk '{
    q=$1; al[q]+=$4; ql[q]=$5
} END{
    for(c in al) {
        printf("%s\t%d\t%d\t%.6f\n", c, al[c], ql[c], al[c]/ql[c])
    }
}' asm_vs_plastid.tsv > alcr_plastid.tsv

awk '{
    q=$1; al[q]+=$4; ql[q]=$5
} END{
    for(c in al) {
        printf("%s\t%d\t%d\t%.6f\n", c, al[c], ql[c], al[c]/ql[c])
    }
}' asm_vs_mito.tsv > alcr_mito.tsv

# -----------------------------
# STEP 4: Classify contigs
# -----------------------------

echo "[STEP 4] Classifying putative organelle-derived contigs"

python3 - << PY
import pandas as pd

ALCR_THRESHOLD = float("${ALCR_THRESHOLD}")
R_NUC_MIN = float("${R_NUC_MIN}")
SDR_PLASTID = float("${SDR_PLASTID}")
SDR_MITO = float("${SDR_MITO}")

depth = pd.read_csv(
    "depth_nuc.tsv",
    sep="\t",
    header=None,
    names=["contig", "D_contig"]
)

alcr_p = pd.read_csv(
    "alcr_plastid.tsv",
    sep="\t",
    header=None,
    names=["contig", "aligned_p", "qlen_p", "ALCR_plastid"]
)

alcr_m = pd.read_csv(
    "alcr_mito.tsv",
    sep="\t",
    header=None,
    names=["contig", "aligned_m", "qlen_m", "ALCR_mito"]
)

df = depth.merge(alcr_p[["contig", "ALCR_plastid"]], on="contig", how="left")
df = df.merge(alcr_m[["contig", "ALCR_mito"]], on="contig", how="left")
df = df.fillna({"ALCR_plastid": 0.0, "ALCR_mito": 0.0})

with open("D_nuc.txt") as f:
    D_nuc = float(f.read().strip())

with open("D_plastid.txt") as f:
    D_plastid = float(f.read().strip())

with open("D_mito.txt") as f:
    D_mito = float(f.read().strip())

rows = []

for _, r in df.iterrows():
    contig = r.contig
    depth_contig = r.D_contig
    alcr_p = r.ALCR_plastid
    alcr_m = r.ALCR_mito

    ratio_to_nuc = depth_contig / D_nuc if D_nuc > 0 else 0.0
    sdr_p = depth_contig / D_plastid if D_plastid > 0 else 0.0
    sdr_m = depth_contig / D_mito if D_mito > 0 else 0.0

    decision = "KEEP"
    organelle_type = "nuclear"
    reason = []

    if alcr_p >= ALCR_THRESHOLD and sdr_p >= SDR_PLASTID and ratio_to_nuc >= R_NUC_MIN:
        decision = "REMOVE_ORGANELLE"
        organelle_type = "plastid"
        reason = ["ALCR_high_plastid", "SDR_high_plastid", "depth_high"]
    elif alcr_m >= ALCR_THRESHOLD and sdr_m >= SDR_MITO and ratio_to_nuc >= R_NUC_MIN:
        decision = "REMOVE_ORGANELLE"
        organelle_type = "mitochondrion"
        reason = ["ALCR_high_mito", "SDR_high_mito", "depth_high"]
    elif (0.3 <= alcr_p < ALCR_THRESHOLD) or (0.3 <= alcr_m < ALCR_THRESHOLD):
        decision = "BORDERLINE"
        reason = ["ALCR_mid"]

    rows.append([
        contig,
        alcr_p,
        alcr_m,
        depth_contig,
        ratio_to_nuc,
        sdr_p,
        sdr_m,
        decision,
        organelle_type,
        ";".join(reason)
    ])

out = pd.DataFrame(
    rows,
    columns=[
        "contig",
        "ALCR_plastid",
        "ALCR_mito",
        "D_contig",
        "ratio_to_nuc",
        "SDR_plastid",
        "SDR_mito",
        "decision",
        "organelle_type",
        "reason"
    ]
)

out.sort_values(["decision", "D_contig"], ascending=[True, False], inplace=True)
out.to_csv("contig_decisions.tsv", sep="\t", index=False)

out[out.decision == "REMOVE_ORGANELLE"]["contig"].to_csv(
    "remove_ids.txt", index=False, header=False
)
out[out.decision == "BORDERLINE"]["contig"].to_csv(
    "borderline_ids.txt", index=False, header=False
)
out[out.decision == "KEEP"]["contig"].to_csv(
    "keep_ids_raw.txt", index=False, header=False
)

print("[INFO] Classification summary")
print(out["decision"].value_counts().to_string())
PY

# -----------------------------
# STEP 5: Produce cleaned assembly
# -----------------------------

echo "[STEP 5] Producing cleaned nuclear assembly"

awk '/^>/{gsub(/^>/,""); split($1,a," "); print a[1]}' "${ASM}" > all_ids_full.txt

if [[ -s remove_ids.txt ]]; then
    grep -v -F -f remove_ids.txt all_ids_full.txt > keep_ids.txt
else
    cp all_ids_full.txt keep_ids.txt
fi

if [[ "${REMOVAL_MODE}" == "remove" ]]; then
    seqtk subseq "${ASM}" keep_ids.txt > nuclear_cleaned.fa
elif [[ "${REMOVAL_MODE}" == "mask" ]]; then
    echo "[FATAL] Masking mode is not implemented in this repository version"
    exit 1
else
    echo "[FATAL] REMOVAL_MODE must be remove or mask"
    exit 1
fi

# -----------------------------
# STEP 6: Summary
# -----------------------------

echo "[STEP 6] Writing summary"

ORIG_N=$(grep -c "^>" "${ASM}")
CLEAN_N=$(grep -c "^>" nuclear_cleaned.fa)
REMOVED_N=$(( ORIG_N - CLEAN_N ))

{
    echo "Organelle filtering summary"
    echo ""
    echo "Original contigs: ${ORIG_N}"
    echo "Cleaned contigs:  ${CLEAN_N}"
    echo "Removed contigs:  ${REMOVED_N}"
    echo ""
    echo "Depth measurements:"
    echo "  Nuclear median: ${D_NUC}x"
    echo "  Plastid mean: ${D_PLASTID}x"
    echo "  Mitochondrial mean: ${D_MITO}x"
    echo ""
    echo "Thresholds used:"
    echo "  Alignment coverage >= ${ALCR_THRESHOLD}"
    echo "  Plastid SDR >= ${SDR_PLASTID}"
    echo "  Mitochondrial SDR >= ${SDR_MITO}"
    echo "  Ratio to nuclear median >= ${R_NUC_MIN}x"
    echo ""
    echo "Key outputs:"
    echo "  contig_decisions.tsv"
    echo "  remove_ids.txt"
    echo "  borderline_ids.txt"
    echo "  nuclear_cleaned.fa"
} | tee summary.txt

echo "[DONE] Organelle filtering completed"