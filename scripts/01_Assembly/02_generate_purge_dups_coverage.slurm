#!/usr/bin/env bash
#SBATCH --job-name=purge_dups_coverage
#SBATCH --output=purge_dups_coverage.%j.out
#SBATCH --error=purge_dups_coverage.%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Generate coverage statistics for purge_dups.
#
# Workflow:
#   1. Map PacBio HiFi reads to the primary assembly using minimap2.
#   2. Calculate read-depth statistics with pbcstat.
#   3. Estimate initial coverage cutoffs using calcuts.
#
# The PB.stat histogram was inspected manually before running
# purge_dups. For the final purge_dups run, the cutoff values
# were manually adjusted to:
#
#   30 80 250

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
READS="/path/to/pacbio_hifi_reads.fastq.gz"
ASSEMBLY="/path/to/hifiasm_primary_assembly.fa"

PREFIX="pb.083025"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# -----------------------------
# Software requirements
# -----------------------------
# purge_dups utilities (pbcstat, calcuts)
# minimap2
# pigz (optional)

for cmd in minimap2 pbcstat calcuts; do
    command -v "${cmd}" >/dev/null || {
        echo "[ERROR] ${cmd} not found."
        exit 127
    }
done

if command -v pigz >/dev/null 2>&1; then
    COMPRESS="pigz -p ${THREADS} -1"
else
    COMPRESS="gzip -1"
fi

echo "[STEP 1] Mapping HiFi reads to the primary assembly"
minimap2 \
    -x map-hifi \
    -t "${THREADS}" \
    "${ASSEMBLY}" \
    "${READS}" \
| ${COMPRESS} \
> "${PREFIX}.paf.gz"

echo "[STEP 2] Running pbcstat"
pbcstat "${PREFIX}.paf.gz"

echo "[STEP 3] Running calcuts"
calcuts PB.stat > cutoffs 2> calcuts.log

# Optional: generate a coverage histogram if hist_plot.py is available.
HIST_PLOT="$(command -v hist_plot.py || true)"
if [[ -n "${HIST_PLOT}" ]]; then
    python "${HIST_PLOT}" PB.stat PB.stat.png || true
fi

echo "[DONE] Coverage files generated: ${PREFIX}.paf.gz, PB.base.cov, PB.stat, cutoffs"