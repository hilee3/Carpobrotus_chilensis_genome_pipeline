#!/usr/bin/env python3

import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate the proportion of well-aligned long reads. "
            "A well-aligned read is defined as having >=90% of its "
            "length aligned in the primary alignment and MAPQ >20."
        )
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output summary file")
    args = parser.parse_args()

    total_primary = 0
    mapped_primary = 0
    well_aligned = 0

    # Additional counts are useful for checking which criterion removes reads.
    aligned_ge90 = 0
    mapq_gt20 = 0

    bam = pysam.AlignmentFile(args.bam, "rb")

    for read in bam.fetch(until_eof=True):

        # Exclude secondary and supplementary alignments.
        # This leaves one primary record per read.
        if read.is_secondary or read.is_supplementary:
            continue

        total_primary += 1

        if read.is_unmapped:
            continue

        mapped_primary += 1

        # Full original read length inferred from the CIGAR.
        # This includes clipped sequence, including hard clipping if present.
        read_length = read.infer_read_length()

        if read_length is None or read_length == 0:
            continue

        # Query bases participating in the primary alignment.
        # Soft-clipped bases are excluded.
        aligned_length = read.query_alignment_length

        aligned_fraction = aligned_length / read_length

        if aligned_fraction >= 0.90:
            aligned_ge90 += 1

        if read.mapping_quality > 20:
            mapq_gt20 += 1

        if aligned_fraction >= 0.90 and read.mapping_quality > 20:
            well_aligned += 1

    bam.close()

    with open(args.output, "w") as out:
        out.write("Metric\tCount\tPercent_of_all_primary_reads\tPercent_of_mapped_primary_reads\n")

        out.write(
            f"Total_primary_reads\t{total_primary}\t100.00\tNA\n"
        )

        out.write(
            f"Mapped_primary_reads\t{mapped_primary}\t"
            f"{100 * mapped_primary / total_primary:.2f}\t100.00\n"
        )

        out.write(
            f"Aligned_fraction_ge_90pct\t{aligned_ge90}\t"
            f"{100 * aligned_ge90 / total_primary:.2f}\t"
            f"{100 * aligned_ge90 / mapped_primary:.2f}\n"
        )

        out.write(
            f"MAPQ_gt_20\t{mapq_gt20}\t"
            f"{100 * mapq_gt20 / total_primary:.2f}\t"
            f"{100 * mapq_gt20 / mapped_primary:.2f}\n"
        )

        out.write(
            f"Well_aligned_ge90pct_and_MAPQ_gt20\t{well_aligned}\t"
            f"{100 * well_aligned / total_primary:.2f}\t"
            f"{100 * well_aligned / mapped_primary:.2f}\n"
        )


if __name__ == "__main__":
    main()
