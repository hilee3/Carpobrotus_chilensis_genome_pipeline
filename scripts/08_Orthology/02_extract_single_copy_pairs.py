import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract single-copy ortholog pairs from OrthoFinder output."
    )
    parser.add_argument(
        "orthogroup_dir",
        type=Path,
        help="Path to the OrthoFinder Orthogroups directory.",
    )
    parser.add_argument(
        "output_tsv",
        type=Path,
        help="Output TSV file.",
    )
    parser.add_argument(
        "--species-a",
        default="chilensis",
        help="Species A column name in Orthogroups.tsv.",
    )
    parser.add_argument(
        "--species-b",
        default="edulis",
        help="Species B column name in Orthogroups.tsv.",
    )
    return parser.parse_args()


args = parse_args()

single_copy_file = args.orthogroup_dir / "Orthogroups_SingleCopyOrthologues.txt"
orthogroups_tsv = args.orthogroup_dir / "Orthogroups.tsv"

single_copy_ogs = set()

with open(single_copy_file) as infile:
    for line in infile:
        og = line.strip()
        if og:
            single_copy_ogs.add(og)

with open(orthogroups_tsv, newline="") as infile, open(
    args.output_tsv, "w", newline=""
) as outfile:

    reader = csv.DictReader(infile, delimiter="\t")

    if args.species_a not in reader.fieldnames or args.species_b not in reader.fieldnames:
        raise SystemExit(
            f"Species columns not found in header: {reader.fieldnames}"
        )

    writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "Orthogroup",
            f"{args.species_a}_gene",
            f"{args.species_b}_gene",
        ]
    )

    n_written = 0

    for row in reader:
        og = row["Orthogroup"]

        if og not in single_copy_ogs:
            continue

        genes_a = [
            x.strip()
            for x in row[args.species_a].split(",")
            if x.strip()
        ]

        genes_b = [
            x.strip()
            for x in row[args.species_b].split(",")
            if x.strip()
        ]

        if len(genes_a) == 1 and len(genes_b) == 1:
            writer.writerow(
                [og, genes_a[0], genes_b[0]]
            )
            n_written += 1

print(f"Saved: {args.output_tsv}")
print(f"Number of pairs: {n_written}")