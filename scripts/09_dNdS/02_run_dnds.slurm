#!/usr/bin/env bash
#SBATCH --job-name=kaks_serial
#SBATCH --output=kaks_serial_%j.out
#SBATCH --error=kaks_serial_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.

WORKDIR="/path/to/genome_comparison"

CHIL_CDS_SOURCE="/path/to/braker.g103t2_removed.longest_isoform.gffread.cds.fa"
EDU_CDS_SOURCE="/path/to/braker.longest_isoform.gffread.cds.fa"

PAIR_TSV="${WORKDIR}/chilensis_edulis_singlecopy_pairs.tsv"

CHIL_PROT="${WORKDIR}/proteomes/chilensis.fa"
EDU_PROT="${WORKDIR}/proteomes/edulis.fa"

CDS_DIR="${WORKDIR}/cds"
PAIR_FASTA_DIR="${WORKDIR}/pair_fastas"
ALN_DIR="${WORKDIR}/alignments"
AXT_DIR="${WORKDIR}/axt"
KAKS_DIR="${WORKDIR}/kaks_results"
SUMMARY_DIR="${WORKDIR}/kaks_summary"

# Directory containing 01_codon_fasta_to_axt.py
SCRIPT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)"

FAIL_LOG="${SUMMARY_DIR}/failed_pairs.tsv"

mkdir -p \
    "${CDS_DIR}" \
    "${PAIR_FASTA_DIR}" \
    "${ALN_DIR}" \
    "${AXT_DIR}" \
    "${KAKS_DIR}" \
    "${SUMMARY_DIR}"

# Copy CDS files into the working directory.
cp \
    "${CHIL_CDS_SOURCE}" \
    "${CDS_DIR}/chilensis_cds.fa"

cp \
    "${EDU_CDS_SOURCE}" \
    "${CDS_DIR}/edulis_cds.fa"

CHIL_CDS="${CDS_DIR}/chilensis_cds.fa"
EDU_CDS="${CDS_DIR}/edulis_cds.fa"

# -----------------------------
# Check required input files
# -----------------------------
for f in \
    "${PAIR_TSV}" \
    "${CHIL_PROT}" \
    "${EDU_PROT}" \
    "${CHIL_CDS}" \
    "${EDU_CDS}" \
    "${SCRIPT_DIR}/01_codon_fasta_to_axt.py"
do
    if [[ ! -s "${f}" ]]; then
        echo "[ERROR] Missing or empty file: ${f}" >&2
        exit 1
    fi
done

# -----------------------------
# Build FASTA indexes
# -----------------------------
rm -f \
    "${CHIL_PROT}.fai" \
    "${EDU_PROT}.fai" \
    "${CHIL_CDS}.fai" \
    "${EDU_CDS}.fai"

samtools faidx "${CHIL_PROT}"
samtools faidx "${EDU_PROT}"
samtools faidx "${CHIL_CDS}"
samtools faidx "${EDU_CDS}"

# Remove CRLF line endings from the pair table if present.
sed -i 's/\r$//' "${PAIR_TSV}"

printf \
    "Orthogroup\tchilensis_gene\tedulis_gene\tfailed_step\n" \
    > "${FAIL_LOG}"

# -----------------------------
# Process ortholog pairs
# -----------------------------
while IFS=$'\t' read -r ORTHOGROUP CHIL_GENE EDU_GENE; do

    # Skip header.
    if [[ "${ORTHOGROUP}" == "Orthogroup" ]]; then
        continue
    fi

    ORTHOGROUP="$(printf '%s' "${ORTHOGROUP}" | tr -d '\r')"
    CHIL_GENE="$(printf '%s' "${CHIL_GENE}" | tr -d '\r')"
    EDU_GENE="$(printf '%s' "${EDU_GENE}" | tr -d '\r')"

    PREFIX="${ORTHOGROUP}"

    PROT_FASTA="${PAIR_FASTA_DIR}/${PREFIX}.aa.fa"
    CDS_FASTA="${PAIR_FASTA_DIR}/${PREFIX}.cds.fa"

    AA_ALN="${ALN_DIR}/${PREFIX}.aa.aln.fa"
    CODON_ALN="${ALN_DIR}/${PREFIX}.codon.aln.fa"

    AXT_FILE="${AXT_DIR}/${PREFIX}.axt"
    KAKS_OUT="${KAKS_DIR}/${PREFIX}.kaks.tsv"

    if [[ -s "${KAKS_OUT}" ]]; then
        continue
    fi

    rm -f \
        "${PROT_FASTA}" \
        "${CDS_FASTA}" \
        "${AA_ALN}" \
        "${CODON_ALN}" \
        "${AXT_FILE}"

    # Fetch C. chilensis protein.
    if ! samtools faidx \
        "${CHIL_PROT}" \
        "${CHIL_GENE}" \
        | sed "1s/^>.*/>chilensis|${CHIL_GENE}/" \
        > "${PROT_FASTA}"
    then
        printf \
            "%s\t%s\t%s\tmissing_chilensis_protein\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${PROT_FASTA}"
        continue
    fi

    # Fetch C. edulis protein.
    if ! samtools faidx \
        "${EDU_PROT}" \
        "${EDU_GENE}" \
        | sed "1s/^>.*/>edulis|${EDU_GENE}/" \
        >> "${PROT_FASTA}"
    then
        printf \
            "%s\t%s\t%s\tmissing_edulis_protein\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${PROT_FASTA}"
        continue
    fi

    # Fetch C. chilensis CDS.
    if ! samtools faidx \
        "${CHIL_CDS}" \
        "${CHIL_GENE}" \
        | sed "1s/^>.*/>chilensis|${CHIL_GENE}/" \
        > "${CDS_FASTA}"
    then
        printf \
            "%s\t%s\t%s\tmissing_chilensis_cds\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${CDS_FASTA}"
        continue
    fi

    # Fetch C. edulis CDS.
    if ! samtools faidx \
        "${EDU_CDS}" \
        "${EDU_GENE}" \
        | sed "1s/^>.*/>edulis|${EDU_GENE}/" \
        >> "${CDS_FASTA}"
    then
        printf \
            "%s\t%s\t%s\tmissing_edulis_cds\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${CDS_FASTA}"
        continue
    fi

    # Protein alignment.
    if ! mafft \
        --quiet \
        "${PROT_FASTA}" \
        > "${AA_ALN}"
    then
        printf \
            "%s\t%s\t%s\tmafft_failed\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${AA_ALN}"
        continue
    fi

    # Codon alignment.
    if ! pal2nal.pl \
        "${AA_ALN}" \
        "${CDS_FASTA}" \
        -output fasta \
        > "${CODON_ALN}"
    then
        printf \
            "%s\t%s\t%s\tpal2nal_failed\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${CODON_ALN}"
        continue
    fi

    # Convert codon FASTA alignment to AXT format.
    if ! python3 \
        "${SCRIPT_DIR}/01_codon_fasta_to_axt.py" \
        "${CODON_ALN}" \
        "${AXT_FILE}"
    then
        printf \
            "%s\t%s\t%s\taxt_conversion_failed\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${AXT_FILE}"
        continue
    fi

    # Calculate Ka, Ks, and Ka/Ks using the YN method.
    if ! KaKs_Calculator \
        -i "${AXT_FILE}" \
        -o "${KAKS_OUT}" \
        -m YN
    then
        printf \
            "%s\t%s\t%s\tkaks_calculator_failed\n" \
            "${ORTHOGROUP}" \
            "${CHIL_GENE}" \
            "${EDU_GENE}" \
            >> "${FAIL_LOG}"

        rm -f "${KAKS_OUT}"
        continue
    fi

done < "${PAIR_TSV}"

# -----------------------------
# Merge pairwise Ka/Ks results
# -----------------------------
SUMMARY_FILE="${SUMMARY_DIR}/all_kaks.tsv"

python3 - << EOF
import glob
import os

kaks_dir = "${KAKS_DIR}"
summary_file = "${SUMMARY_FILE}"

rows = []

for fp in sorted(
    glob.glob(
        os.path.join(
            kaks_dir,
            "*.kaks.tsv"
        )
    )
):
    og = (
        os.path.basename(fp)
        .replace(".kaks.tsv", "")
    )

    with open(fp) as f:
        lines = [
            line.strip()
            for line in f
            if line.strip()
        ]

    if len(lines) < 2:
        continue

    fields = lines[1].split("\t")

    if len(fields) < 5:
        continue

    rows.append(
        [og] + fields
    )

with open(summary_file, "w") as out:
    out.write(
        "Orthogroup\t"
        "Sequence\t"
        "Method\t"
        "Ka\t"
        "Ks\t"
        "Ka_Ks\t"
        "P_Value(Fisher)\t"
        "Length\t"
        "S-Sites\t"
        "N-Sites\t"
        "Fold-Sites(0:2:4)\t"
        "Substitutions\t"
        "S-Substitutions\t"
        "N-Substitutions\t"
        "Fold-S-Substitutions\t"
        "Fold-N-Substitutions\t"
        "Divergence-Time\t"
        "Substitution-Rate-Ratio"
        "(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\n"
    )

    for row in rows:
        out.write(
            "\t".join(row) + "\n"
        )

print(f"[INFO] Summary written to: {summary_file}")
print(f"[INFO] Total valid results: {len(rows)}")
EOF

echo "[INFO] Ka/Ks analysis completed."
echo "[INFO] Summary: ${SUMMARY_FILE}"
echo "[INFO] Failure log: ${FAIL_LOG}"