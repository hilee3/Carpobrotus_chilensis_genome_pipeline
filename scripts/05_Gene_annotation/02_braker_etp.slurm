#!/bin/bash
#SBATCH --job-name=braker_chilensis_etp
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --output=braker_%j.out
#SBATCH --error=braker_%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
GENOME="/path/to/softmasked_genome.fa"
BAM="/path/to/rnaseq.sorted.bam"
PROT_SEQ="/path/to/protein_database.fa"
FINAL_OUT="/path/to/output/braker3"

# Keep species name stable unless you intentionally want a new training set every run
SPECIES="carpobrotus_chilensis_etp_contig178"

# -----------------------------
# Scratch layout
# -----------------------------
BASE_SCRATCH="${SCRATCH:-/scratch/$USER}"
JOBTMP="${BASE_SCRATCH}/braker_${SLURM_JOB_ID:-$$}"
WORK="${JOBTMP}/work"
ARCHIVE_NAME="braker_workdir_${SLURM_JOB_ID:-$$}.tar.gz"

mkdir -p "$WORK" "$FINAL_OUT"

export TMPDIR="$JOBTMP"
export TEMP="$JOBTMP"
export TMP="$JOBTMP"

# -----------------------------
# Software requirements
# -----------------------------
# Activate an environment containing:
#   BRAKER3
#   AUGUSTUS
#   GeneMark-ETP
#   ProtHint
#   DIAMOND
#   StringTie
#   samtools
#   bamtools
#   bedtools
#   gffread
#
# GeneMark requires a valid license key, usually at:
#   ~/.gm_key

# GeneMark-ETP installation
GMETP_BASE="/path/to/gmetp_linux_64"
export GENEMARK_PATH="${GMETP_BASE}/bin"
export PATH="${GMETP_BASE}/bin:${GMETP_BASE}/bin/gmes:${GMETP_BASE}/tools:${PATH}"
export GM_KEY="${HOME}/.gm_key"

# ProHint installation
PROTHINT_BASE="/path/to/ProtHint"
export PROTHINT_PATH="${PROTHINT_BASE}/bin"
export PATH="${PROTHINT_BASE}/bin:${PATH}"

# -----------------------------
# Check inputs and dependencies
# -----------------------------
[[ -r "$GENOME" ]] || { echo "[ERR] Genome not found: $GENOME"; exit 1; }
[[ -r "$BAM" ]] || { echo "[ERR] BAM not found: $BAM"; exit 1; }
[[ -r "${BAM}.bai" || -r "${BAM%.bam}.bai" ]] || { echo "[ERR] BAM index not found"; exit 1; }
[[ -r "$PROT_SEQ" ]] || { echo "[ERR] Protein DB not found: $PROT_SEQ"; exit 1; }
[[ -f "$GM_KEY" ]] || { echo "[ERR] GeneMark license not found at $GM_KEY"; exit 1; }

command -v braker.pl >/dev/null || { echo "[ERR] braker.pl not in PATH"; exit 1; }
command -v augustus >/dev/null || { echo "[ERR] augustus not in PATH"; exit 1; }
command -v samtools >/dev/null || { echo "[ERR] samtools not in PATH"; exit 1; }
command -v bamtools >/dev/null || { echo "[ERR] bamtools not in PATH"; exit 1; }
command -v diamond >/dev/null || { echo "[ERR] diamond not in PATH"; exit 1; }
command -v stringtie >/dev/null || { echo "[ERR] stringtie not in PATH"; exit 1; }
command -v bedtools >/dev/null || { echo "[ERR] bedtools not in PATH"; exit 1; }
command -v gffread >/dev/null || { echo "[ERR] gffread not in PATH"; exit 1; }
command -v rsync >/dev/null || { echo "[ERR] rsync not in PATH"; exit 1; }

# GeneMark executables
command -v gmetp.pl >/dev/null || { echo "[ERR] gmetp.pl not in PATH"; exit 1; }
command -v gmes_petap.pl >/dev/null || { echo "[ERR] gmes_petap.pl not in PATH"; exit 1; }

# ProHint executables
command -v prothint.py >/dev/null || { echo "[ERR] prothint.py not in PATH"; exit 1; }

# -----------------------------
# Writable AUGUSTUS config
# -----------------------------
AUG_BIN="$(command -v augustus || true)"
if [[ -n "${AUG_BIN}" && -z "${AUGUSTUS_CONFIG_PATH:-}" ]]; then
    export AUGUSTUS_CONFIG_PATH="$(dirname "$AUG_BIN")/../config"
fi

if [[ -z "${AUGUSTUS_CONFIG_PATH:-}" ]]; then
    echo "[ERR] AUGUSTUS_CONFIG_PATH could not be determined"
    exit 1
fi

if [[ ! -d "${AUGUSTUS_CONFIG_PATH}/species" ]]; then
    echo "[ERR] AUGUSTUS_CONFIG_PATH does not look valid: ${AUGUSTUS_CONFIG_PATH}"
    exit 1
fi

if [[ ! -w "${AUGUSTUS_CONFIG_PATH}/species" ]]; then
    echo "[INFO] AUGUSTUS config is not writable; copying to scratch"
    LOCAL_CFG="${JOBTMP}/augustus_config"
    rsync -a "${AUGUSTUS_CONFIG_PATH}/" "${LOCAL_CFG}/"
    export AUGUSTUS_CONFIG_PATH="$LOCAL_CFG"
fi

# -----------------------------
# Cleanup / copy-back
# -----------------------------
finish() {
    set +e

    echo "[INFO] Syncing BRAKER outputs to ${FINAL_OUT}"

    [[ -e "${WORK}/braker.gtf" ]]        && rsync -a "${WORK}/braker.gtf" "${FINAL_OUT}/"
    [[ -e "${WORK}/braker.gff3" ]]       && rsync -a "${WORK}/braker.gff3" "${FINAL_OUT}/"
    [[ -e "${WORK}/braker.aa" ]]         && rsync -a "${WORK}/braker.aa" "${FINAL_OUT}/"
    [[ -e "${WORK}/braker.codingseq" ]]  && rsync -a "${WORK}/braker.codingseq" "${FINAL_OUT}/"
    [[ -e "${WORK}/hintsfile.gff" ]]     && rsync -a "${WORK}/hintsfile.gff" "${FINAL_OUT}/"
    [[ -e "${WORK}/what-to-cite.txt" ]]  && rsync -a "${WORK}/what-to-cite.txt" "${FINAL_OUT}/"
    [[ -e "${WORK}/braker.log" ]]        && rsync -a "${WORK}/braker.log" "${FINAL_OUT}/"

    [[ -d "${WORK}/GeneMark-ETP" ]]      && rsync -a "${WORK}/GeneMark-ETP" "${FINAL_OUT}/"
    [[ -d "${WORK}/Augustus" ]]          && rsync -a "${WORK}/Augustus" "${FINAL_OUT}/"
    [[ -d "${WORK}/errors" ]]            && rsync -a "${WORK}/errors" "${FINAL_OUT}/"
    [[ -d "${WORK}/bbc" ]]               && rsync -a "${WORK}/bbc" "${FINAL_OUT}/"
    [[ -d "${WORK}/braker_original" ]]   && rsync -a "${WORK}/braker_original" "${FINAL_OUT}/"

    echo "[INFO] Archiving work directory"
    tar -czf "${FINAL_OUT}/${ARCHIVE_NAME}" -C "${BASE_SCRATCH}" "$(basename "${JOBTMP}")"

    echo "[INFO] Cleaning scratch"
    rm -rf "${JOBTMP}"
}
trap finish EXIT INT TERM

# -----------------------------
# Run BRAKER ETP
# -----------------------------
echo "[INFO] Host: $(hostname)"
echo "[INFO] Work dir: ${WORK}"
echo "[INFO] Final out: ${FINAL_OUT}"
echo "[INFO] AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH}"
echo "[INFO] GENEMARK_PATH=${GENEMARK_PATH}"

braker.pl \
    --genome="$GENOME" \
    --prot_seq="$PROT_SEQ" \
    --bam="$BAM" \
    --species="$SPECIES" \
    --softmasking \
    --threads="${SLURM_CPUS_PER_TASK:-16}" \
    --workingdir="$WORK" \
    --GENEMARK_PATH="$GENEMARK_PATH" \
    --gff3
