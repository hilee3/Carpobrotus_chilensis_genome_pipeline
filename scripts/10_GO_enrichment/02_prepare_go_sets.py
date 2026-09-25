import csv
import math
from pathlib import Path


# -----------------------------
# Inputs
# -----------------------------

WORKDIR = Path(
    "/path/to/genome_comparison"
)

KAKS_FILTERED = (
    WORKDIR
    / "kaks_summary"
    / "all_kaks.filtered.tsv"
)

OUTDIR = (
    WORKDIR
    / "go_enrichment"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


study_out = (
    OUTDIR
    / "dnds_gt1_ds001_genes.txt"
)

background_out = (
    OUTDIR
    / "kaks_background_genes.txt"
)

summary_out = (
    OUTDIR
    / "dnds_gt1_ds001_selected_rows.tsv"
)


def clean_gene_id(
    seq_id,
):
    """
    Convert an identifier such as
    chilensis|g11487.t1 to g11487.t1.
    """

    return (
        seq_id
        .strip()
        .split("|")[-1]
    )


background = []
study = []
study_rows = []


with open(
    KAKS_FILTERED,
    newline="",
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    required_cols = {
        "Sequence",
        "Ka_Ks",
        "Ks",
    }

    missing = (
        required_cols
        - set(
            reader.fieldnames
            or []
        )
    )

    if missing:
        raise SystemExit(
            f"Missing required columns: {missing}"
        )

    for row in reader:

        raw_gene = row[
            "Sequence"
        ]

        gene = clean_gene_id(
            raw_gene
        )

        if not gene:
            continue

        background.append(
            gene
        )

        try:
            dnds = float(
                row["Ka_Ks"]
            )

            ds = float(
                row["Ks"]
            )

        except ValueError:
            continue

        if not (
            math.isfinite(dnds)
            and math.isfinite(ds)
        ):
            continue

        # Study set:
        # dN/dS > 1 and dS >= 0.01.
        if (
            dnds > 1
            and ds >= 0.01
        ):
            study.append(
                gene
            )

            row[
                "clean_gene_id"
            ] = gene

            study_rows.append(
                row
            )


# Remove duplicates while preserving order.
background = list(
    dict.fromkeys(
        background
    )
)

study = list(
    dict.fromkeys(
        study
    )
)


with open(
    background_out,
    "w",
) as outfile:

    for gene in background:
        outfile.write(
            gene + "\n"
        )


with open(
    study_out,
    "w",
) as outfile:

    for gene in study:
        outfile.write(
            gene + "\n"
        )


if study_rows:

    fieldnames = list(
        study_rows[0].keys()
    )

    with open(
        summary_out,
        "w",
        newline="",
    ) as outfile:

        writer = csv.DictWriter(
            outfile,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        writer.writerows(
            study_rows
        )

print(f"Input file: {KAKS_FILTERED}")
print(f"Background genes: {len(background)}")
print(f"Study genes: {len(study)}")
print(f"Background file: {background_out}")
print(f"Study file: {study_out}")
print(f"Selected rows: {summary_out}")