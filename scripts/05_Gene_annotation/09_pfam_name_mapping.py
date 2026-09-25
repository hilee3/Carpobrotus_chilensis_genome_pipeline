#!/usr/bin/env python3

import csv
import json
import time
import urllib.request
import urllib.error


# -----------------------------
# Inputs
# -----------------------------

# Input table containing the historical Pfam domains reported by
# Piriyapongsa et al. (2007). Required columns: "Accession" and "ID".
input_file = "/path/to/Piriyapongsa_2007_124_TE_associated_Pfam_domains.tsv"
output_file = "/path/to/Piriyapongsa_Pfam_name_mapping.tsv"


# -----------------------------
# Map historical Pfam names to
# current names using accessions
# -----------------------------

rows = []

with open(input_file) as f:
    reader = csv.DictReader(
        f,
        delimiter="\t",
    )

    for row in reader:
        accession = row["Accession"]
        old_name = row["ID"]

        url = (
            "https://www.ebi.ac.uk/interpro/api/"
            f"entry/pfam/{accession}/"
        )

        current_name = ""
        status = ""

        try:
            with urllib.request.urlopen(
                url,
                timeout=30,
            ) as response:
                data = json.load(response)

            current_name = (
                data["metadata"]["name"]["short"]
            )

            status = "resolved"

        except urllib.error.HTTPError as e:
            status = f"HTTP_{e.code}"

        except Exception as e:
            status = f"ERROR:{type(e).__name__}"

        rows.append(
            [
                accession,
                old_name,
                current_name,
                status,
            ]
        )

        print(
            accession,
            old_name,
            current_name,
            status,
        )

        time.sleep(0.2)


# -----------------------------
# Write mapping table
# -----------------------------

with open(
    output_file,
    "w",
    newline="",
) as f:

    writer = csv.writer(
        f,
        delimiter="\t",
    )

    writer.writerow(
        [
            "Accession",
            "Old_Pfam_name",
            "Current_Pfam_name",
            "Status",
        ]
    )

    writer.writerows(rows)

print(f"Wrote: {output_file}")
