#!/usr/bin/env bash
#SBATCH --job-name=edta_edulis
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --output=edta_edulis_%j.out
#SBATCH --error=edta_edulis_%j.err

set -euo pipefail
umask 002
ulimit -n 65535

# Required software:
#   EDTA
#   RepeatMasker

# -----------------------------
# Inputs
# -----------------------------

GENOME_SRC="/path/to/unmasked_nuclear_only.fa"
WORKDIR="/path/to/output/edulis_edta"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
SPECIES="others"

# Scratch location
SCRATCH_BASE="/path/to/scratch"
RUN="${SCRATCH_BASE}/edta_${SLURM_JOB_ID:-$$}"
WORK="${RUN}/run"
TMP="${RUN}/tmp"

mkdir -p "${WORK}" "${TMP}" "${WORKDIR}"
export TMPDIR="${TMP}"

# EDTA writes output files next to the input genome.
cp "${GENOME_SRC}" "${WORK}/genome.fa"

cd "${WORK}"

# -----------------------------
# EDTA
# -----------------------------

EDTA.pl \
    --genome genome.fa \
    --species "${SPECIES}" \
    --step all \
    --anno 1 \
    --threads "${THREADS}" \
    > EDTA.stdout.log \
    2> EDTA.stderr.log

# -----------------------------
# RepeatMasker soft masking
# using the EDTA TE library
# -----------------------------

TELIB="$(ls -1 *EDTA.TElib.fa 2>/dev/null | head -n 1 || true)"

if [[ -z "${TELIB}" || ! -s "${TELIB}" ]]; then
    echo "[ERR] EDTA TE library not found"
    exit 1
fi

RepeatMasker \
    -pa "${THREADS}" \
    -xsmall \
    -lib "${TELIB}" \
    -dir . \
    -gff \
    genome.fa \
    > RepeatMasker.stdout.log \
    2> RepeatMasker.stderr.log

# -----------------------------
# Export outputs
# -----------------------------

for pat in \
    "*EDTA.TElib.fa" \
    "*EDTA.TElib.novel.fa" \
    "*EDTA.intact.gff3" \
    "*EDTA.TEanno.gff3" \
    "*EDTA.TEanno.sum" \
    "*TEanno.density_plots.pdf" \
    "*divergence_plot.pdf" \
    "genome.fa.masked" \
    "genome.fa.out" \
    "genome.fa.tbl" \
    "genome.fa.out.gff" \
    "EDTA.log" \
    "EDTA.stdout.log" \
    "EDTA.stderr.log" \
    "RepeatMasker.stdout.log" \
    "RepeatMasker.stderr.log"
do
    for f in ${pat}; do
        [[ -e "${f}" ]] && cp -f "${f}" "${WORKDIR}/"
    done
done

if [[ -s genome.fa.masked ]]; then
    cp -f genome.fa.masked \
        "${WORKDIR}/edulis_nuclear_only.softmasked.fa"
fi

echo "[INFO] Output: ${WORKDIR}"
