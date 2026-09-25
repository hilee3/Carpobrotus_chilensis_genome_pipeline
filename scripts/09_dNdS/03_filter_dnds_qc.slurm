#!/bin/bash
#SBATCH --job-name=filter_kaks_qc
#SBATCH --output=filter_kaks_qc_%j.out
#SBATCH --error=filter_kaks_qc_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/genome_comparison"

export WORKDIR

python3 - << 'PYEOF'
import csv
import math
import os
from pathlib import Path
from statistics import mean, median

from Bio import SeqIO


workdir = Path(os.environ["WORKDIR"])

input_summary = (
    workdir
    / "kaks_summary"
    / "all_kaks.tsv"
)

align_dir = (
    workdir
    / "alignments"
)

output_filtered = (
    workdir
    / "kaks_summary"
    / "all_kaks.filtered.tsv"
)

output_removed = (
    workdir
    / "kaks_summary"
    / "filtered_out.tsv"
)

output_stats = (
    workdir
    / "kaks_summary"
    / "filter_stats.txt"
)


def find_column(
    fieldnames,
    candidates,
):
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate

    return None


def ungapped_length(seq):
    return len(
        str(seq)
        .replace("-", "")
        .replace(".", "")
    )


def has_internal_stop_codon(
    codon_fasta,
):
    """
    Return True if any sequence contains an internal stop codon.

    A terminal stop codon is allowed.
    """

    records = list(
        SeqIO.parse(
            str(codon_fasta),
            "fasta",
        )
    )

    if len(records) == 0:
        return True

    stop_codons = {
        "TAA",
        "TAG",
        "TGA",
    }

    for record in records:
        seq = (
            str(record.seq)
            .upper()
            .replace("-", "")
            .replace(".", "")
        )

        if len(seq) < 3:
            return True

        if len(seq) % 3 != 0:
            return True

        codons = [
            seq[i:i + 3]
            for i in range(
                0,
                len(seq),
                3,
            )
        ]

        for codon in codons[:-1]:
            if codon in stop_codons:
                return True

    return False


with open(
    input_summary,
    "r",
    newline="",
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    fieldnames = reader.fieldnames

    if fieldnames is None:
        raise SystemExit(
            "Input summary file has no header."
        )

    ka_col = find_column(
        fieldnames,
        ["Ka", "dN"],
    )

    ks_col = find_column(
        fieldnames,
        ["Ks", "dS"],
    )

    kaks_col = find_column(
        fieldnames,
        [
            "Ka_Ks",
            "Ka/Ks",
            "dN_dS",
            "omega",
        ],
    )

    og_col = find_column(
        fieldnames,
        [
            "Orthogroup",
            "orthogroup",
        ],
    )

    if (
        ka_col is None
        or ks_col is None
        or kaks_col is None
        or og_col is None
    ):
        raise SystemExit(
            "Could not identify required columns. "
            f"Found header: {fieldnames}"
        )

    kept_rows = []
    removed_rows = []

    total = 0

    removed_non_numeric = 0
    removed_ds_range = 0
    removed_kaks_gt_10 = 0
    removed_missing_alignment = 0
    removed_short_aa = 0
    removed_short_codon = 0
    removed_internal_stop = 0

    kept_ka = []
    kept_ks = []
    kept_kaks = []

    for row in reader:
        total += 1

        og = row[og_col]

        try:
            ka = float(
                row[ka_col]
            )

            ks = float(
                row[ks_col]
            )

            kaks = float(
                row[kaks_col]
            )

        except (
            TypeError,
            ValueError,
        ):
            row[
                "filter_reason"
            ] = "non_numeric_value"

            removed_rows.append(
                row
            )

            removed_non_numeric += 1
            continue

        if not (
            math.isfinite(ka)
            and math.isfinite(ks)
            and math.isfinite(kaks)
        ):
            row[
                "filter_reason"
            ] = "non_finite_value"

            removed_rows.append(
                row
            )

            removed_non_numeric += 1
            continue

        # Retain 0.001 <= dS <= 2.
        if ks < 0.001 or ks > 2:
            row[
                "filter_reason"
            ] = "dS_out_of_range"

            removed_rows.append(
                row
            )

            removed_ds_range += 1
            continue

        # Retain dN/dS <= 10.
        if kaks > 10:
            row[
                "filter_reason"
            ] = "dN_dS_gt_10"

            removed_rows.append(
                row
            )

            removed_kaks_gt_10 += 1
            continue

        aa_aln = (
            align_dir
            / f"{og}.aa.aln.fa"
        )

        codon_aln = (
            align_dir
            / f"{og}.codon.aln.fa"
        )

        if (
            not aa_aln.exists()
            or not codon_aln.exists()
        ):
            row[
                "filter_reason"
            ] = "missing_alignment_file"

            removed_rows.append(
                row
            )

            removed_missing_alignment += 1
            continue

        aa_records = list(
            SeqIO.parse(
                str(aa_aln),
                "fasta",
            )
        )

        codon_records = list(
            SeqIO.parse(
                str(codon_aln),
                "fasta",
            )
        )

        if (
            len(aa_records) < 2
            or len(codon_records) < 2
        ):
            row[
                "filter_reason"
            ] = "invalid_alignment_record_count"

            removed_rows.append(
                row
            )

            removed_missing_alignment += 1
            continue

        # Use the shorter ungapped sequence in each pair.
        aa_len = min(
            ungapped_length(
                record.seq
            )
            for record in aa_records
        )

        codon_len = min(
            ungapped_length(
                record.seq
            )
            for record in codon_records
        )

        if aa_len < 50:
            row[
                "filter_reason"
            ] = "aa_alignment_lt_50aa"

            removed_rows.append(
                row
            )

            removed_short_aa += 1
            continue

        if codon_len < 150:
            row[
                "filter_reason"
            ] = "codon_alignment_lt_150nt"

            removed_rows.append(
                row
            )

            removed_short_codon += 1
            continue

        if has_internal_stop_codon(
            codon_aln
        ):
            row[
                "filter_reason"
            ] = "internal_stop_codon"

            removed_rows.append(
                row
            )

            removed_internal_stop += 1
            continue

        row[
            "aa_alignment_length"
        ] = aa_len

        row[
            "codon_alignment_length"
        ] = codon_len

        kept_rows.append(
            row
        )

        kept_ka.append(
            ka
        )

        kept_ks.append(
            ks
        )

        kept_kaks.append(
            kaks
        )


filtered_fieldnames = list(
    fieldnames
)

for column in [
    "aa_alignment_length",
    "codon_alignment_length",
]:
    if column not in filtered_fieldnames:
        filtered_fieldnames.append(
            column
        )

with open(
    output_filtered,
    "w",
    newline="",
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=filtered_fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(
        kept_rows
    )


removed_fieldnames = list(
    fieldnames
)

if (
    "filter_reason"
    not in removed_fieldnames
):
    removed_fieldnames.append(
        "filter_reason"
    )

with open(
    output_removed,
    "w",
    newline="",
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=removed_fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(
        removed_rows
    )


with open(
    output_stats,
    "w",
) as outfile:

    outfile.write(
        f"Input summary\t{input_summary}\n"
    )

    outfile.write(
        f"Alignment directory\t{align_dir}\n"
    )

    outfile.write(
        f"Filtered output\t{output_filtered}\n"
    )

    outfile.write(
        f"Removed output\t{output_removed}\n"
    )

    outfile.write(
        f"Total rows\t{total}\n"
    )

    outfile.write(
        f"Kept rows\t{len(kept_rows)}\n"
    )

    outfile.write(
        f"Removed rows\t{len(removed_rows)}\n"
    )

    outfile.write("\n")

    outfile.write(
        "Applied filters\n"
    )

    outfile.write(
        "0.001 <= dS <= 2\n"
    )

    outfile.write(
        "dN/dS <= 10\n"
    )

    outfile.write(
        "codon alignment length >= 150 nt\n"
    )

    outfile.write(
        "amino acid alignment length >= 50 aa\n"
    )

    outfile.write(
        "remove alignments containing internal stop codons\n"
    )

    outfile.write("\n")

    outfile.write(
        "Removed non-numeric or non-finite\t"
        f"{removed_non_numeric}\n"
    )

    outfile.write(
        "Removed dS out of range\t"
        f"{removed_ds_range}\n"
    )

    outfile.write(
        "Removed dN/dS > 10\t"
        f"{removed_kaks_gt_10}\n"
    )

    outfile.write(
        "Removed missing/invalid alignment\t"
        f"{removed_missing_alignment}\n"
    )

    outfile.write(
        "Removed AA alignment < 50 aa\t"
        f"{removed_short_aa}\n"
    )

    outfile.write(
        "Removed codon alignment < 150 nt\t"
        f"{removed_short_codon}\n"
    )

    outfile.write(
        "Removed internal stop codon\t"
        f"{removed_internal_stop}\n"
    )

    if kept_rows:
        outfile.write("\n")

        outfile.write(
            f"Mean Ka\t{mean(kept_ka):.6f}\n"
        )

        outfile.write(
            f"Median Ka\t{median(kept_ka):.6f}\n"
        )

        outfile.write(
            f"Mean dS\t{mean(kept_ks):.6f}\n"
        )

        outfile.write(
            f"Median dS\t{median(kept_ks):.6f}\n"
        )

        outfile.write(
            f"Mean dN/dS\t{mean(kept_kaks):.6f}\n"
        )

        outfile.write(
            f"Median dN/dS\t{median(kept_kaks):.6f}\n"
        )
PYEOF
