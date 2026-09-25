#!/usr/bin/env bash
#SBATCH --job-name=plot_kaks
#SBATCH --output=plot_kaks_%j.out
#SBATCH --error=plot_kaks_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

WORKDIR="/path/to/genome_comparison"
INPUT_TSV="${WORKDIR}/kaks_summary/all_kaks.filtered.tsv"
OUTDIR="${WORKDIR}/kaks_summary/plots"

mkdir -p "${OUTDIR}"

export INPUT_TSV
export OUTDIR

python3 - << 'PYEOF'
import os

import matplotlib.pyplot as plt
import pandas as pd


input_tsv = os.environ["INPUT_TSV"]
outdir = os.environ["OUTDIR"]


# -----------------------------
# Read QC-filtered Ka/Ks results
# -----------------------------
df = pd.read_csv(
    input_tsv,
    sep="\t",
)

df["Ks"] = pd.to_numeric(
    df["Ks"],
    errors="coerce",
)

df["Ka_Ks"] = pd.to_numeric(
    df["Ka_Ks"],
    errors="coerce",
)

ds = (
    df["Ks"]
    .replace([float("inf"), float("-inf")], pd.NA)
    .dropna()
)

dnds = (
    df["Ka_Ks"]
    .replace([float("inf"), float("-inf")], pd.NA)
    .dropna()
)


# -----------------------------
# dS distribution
# Publication-style version
# -----------------------------
median_ds = ds.median()

fig, ax = plt.subplots(
    figsize=(5.2, 4.2)
)

ax.hist(
    ds,
    bins=100,
    range=(0, 0.5),
    color="steelblue",
    edgecolor="none",
)

ax.axvline(
    median_ds,
    linestyle="--",
    linewidth=1.0,
    color="black",
)

ax.text(
    median_ds + 0.01,
    ax.get_ylim()[1] * 0.7,
    f"Median = {median_ds:.3f}",
    ha="left",
    va="top",
    fontsize=10,
)

ax.set_xlabel(
    "dS",
    fontsize=11,
)

ax.set_ylabel(
    "Number of orthologous gene pairs",
    fontsize=11,
)

ax.set_xlim(
    0,
    0.5,
)

ax.set_xticks(
    [
        0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
    ]
)

ax.tick_params(
    axis="both",
    labelsize=10,
)

plt.tight_layout()

plt.savefig(
    f"{outdir}/dS_histogram.png",
    dpi=600,
    bbox_inches="tight",
)

plt.savefig(
    f"{outdir}/dS_histogram.pdf",
    bbox_inches="tight",
)

plt.close()


# -----------------------------
# dN/dS distribution
# -----------------------------
fig, ax = plt.subplots(
    figsize=(6, 4.5)
)

ax.hist(
    dnds,
    bins=60,
)

ax.axvline(
    x=1,
    linestyle="--",
    linewidth=1.0,
    color="black",
)

ax.set_xlabel(
    "dN/dS"
)

ax.set_ylabel(
    "Count"
)

ax.set_title(
    "Distribution of dN/dS ratios"
)

plt.tight_layout()

plt.savefig(
    f"{outdir}/dNdS_histogram.png",
    dpi=300,
)

plt.close()

PYEOF

echo "[INFO] Plots written to: ${OUTDIR}"