#!/bin/bash
#SBATCH --job-name=edta_repeat_annotation
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --output=edta_repeat_annotation.%j.out
#SBATCH --error=edta_repeat_annotation.%j.err

set -euo pipefail
ulimit -n 65535

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
WORKDIR="/path/to/output/edta"
GENOME_SRC="/path/to/assembly_split.fa"
CDS_FASTA="/path/to/braker.longest_isoform.gffread.cds.fa"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
SPECIES="others"

# Scratch directory
SCRATCH_BASE="${SCRATCH:-/scratch/$USER}"

# -----------------------------
# Software requirements
# -----------------------------
# EDTA (v2.2.2)
# RepeatMasker

# -----------------------------
# Scratch workspace
# -----------------------------
RUN="${SCRATCH_BASE}/edta_${SLURM_JOB_ID:-manual}"
WORK="${RUN}/run"
TMP="${RUN}/tmp"

mkdir -p "$WORK" "$TMP" "$WORKDIR"
export TMPDIR="$TMP"

# Put input genome in scratch because EDTA writes outputs
# next to the genome file.
cp "$GENOME_SRC" "$WORK/genome.fa"

# Save a copy of this script for reproducibility.
cp "$0" "$WORKDIR/$(basename "$0")"

echo "[INFO] Starting EDTA in scratch: $WORK"
echo "[INFO] Genome source: $GENOME_SRC"
echo "[INFO] CDS source: $CDS_FASTA"
echo "[INFO] Threads: $THREADS"

# -----------------------------
# EDTA run
# -----------------------------
cd "$WORK"

EDTA.pl \
  --genome genome.fa \
  --species "$SPECIES" \
  --step all \
  --cds "$CDS_FASTA" \
  --anno 1 \
  --threads "$THREADS" \
  > EDTA.stdout.log \
  2> EDTA.stderr.log

# -----------------------------
# RepeatMasker softmask with EDTA library
# -----------------------------
TELIB="$(ls -1 *EDTA.TElib.fa 2>/dev/null | head -n 1 || true)"

if [[ -n "${TELIB:-}" && -s "$TELIB" ]]; then
  echo "[INFO] Running RepeatMasker -xsmall with EDTA TElib: $TELIB"

  RepeatMasker \
    -pa "$THREADS" \
    -xsmall \
    -lib "$TELIB" \
    -dir . \
    -gff \
    genome.fa \
    > RepeatMasker.stdout.log \
    2> RepeatMasker.stderr.log
else
  echo "[WARN] No EDTA TElib found. Skipping RepeatMasker softmask step." \
    | tee RepeatMasker.stderr.log
fi

# -----------------------------
# Export selected artifacts
# -----------------------------
echo "[INFO] Exporting selected artifacts to $WORKDIR"

BASENAME="$(basename "$GENOME_SRC")"
ROOT="${BASENAME%.*}"

ARTIFACTS=(
  "*EDTA.TElib.fa"
  "*EDTA.TElib.novel.fa"
  "*EDTA.intact.gff3"
  "*EDTA.TEanno.gff3"
  "*EDTA.TEanno.sum"
  "*TEanno.density_plots.pdf"
  "*divergence_plot.pdf"
  "genome.fa.masked"
  "genome.fa.out"
  "genome.fa.tbl"
  "genome.fa.out.gff"
  "EDTA.log"
  "EDTA.stdout.log"
  "EDTA.stderr.log"
  "RepeatMasker.stdout.log"
  "RepeatMasker.stderr.log"
)

for pat in "${ARTIFACTS[@]}"; do
  for f in $pat; do
    [[ -e "$f" ]] && cp -f "$f" "$WORKDIR"/
  done
done

# Friendlier softmasked filename
if [[ -s "genome.fa.masked" ]]; then
  cp -f "genome.fa.masked" "${WORKDIR}/${ROOT}.softmasked.fa"
fi

# -----------------------------
# Scratch cleanup
# -----------------------------
echo "[INFO] Cleaning scratch: $RUN"
rm -rf "$RUN"

echo "[INFO] Done."