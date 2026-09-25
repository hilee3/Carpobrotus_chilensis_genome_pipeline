#!/bin/bash

set -euo pipefail

# Apply corrections for two internal adaptor-like regions
# identified by NCBI FCS-GX in ptg000006l_1:
#   6319403-6319427
#   6340302-6340326

ASSEMBLY="/path/to/cp_cov_filter.filtered_nuclear.tiara_filtered.fa"

# Extract the three retained regions.
samtools faidx "${ASSEMBLY}" \
    ptg000006l_1:1-6319402 \
    > ptg000006l_1_part1.fa

samtools faidx "${ASSEMBLY}" \
    ptg000006l_1:6319428-6340301 \
    > ptg000006l_1_part2.fa

samtools faidx "${ASSEMBLY}" \
    ptg000006l_1:6340327-88555222 \
    > ptg000006l_1_part3.fa

# Rename split-contig headers.
sed -i '1s|.*|>ptg000006l_1a|' ptg000006l_1_part1.fa
sed -i '1s|.*|>ptg000006l_1b|' ptg000006l_1_part2.fa
sed -i '1s|.*|>ptg000006l_1c|' ptg000006l_1_part3.fa

# Remove the original contig.
seqkit grep -v -p "ptg000006l_1" \
    "${ASSEMBLY}" \
    > assembly_without_ptg000006l_1.fa

# Add the three corrected fragments.
cat \
    assembly_without_ptg000006l_1.fa \
    ptg000006l_1_part1.fa \
    ptg000006l_1_part2.fa \
    ptg000006l_1_part3.fa \
    > assembly_split.fa

echo "Final contig count: $(grep -c '^>' assembly_split.fa)"
echo "Done: assembly_split.fa"
