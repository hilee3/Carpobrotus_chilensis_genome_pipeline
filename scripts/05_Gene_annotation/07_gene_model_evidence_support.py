#!/usr/bin/env python3

gtf_file = "/path/to/braker.gtf"
hints_file = "/path/to/hintsfile.gff"

rna_introns = set()
protein_introns = set()

with open(hints_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        if len(fields) < 9 or fields[2] != "intron":
            continue

        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]

        intron = (chrom, start, end, strand)

        if "src=E" in attrs:
            rna_introns.add(intron)

        if (
            "src=P" in attrs or
            "src=C" in attrs or
            "src=M" in attrs
        ):
            protein_introns.add(intron)


gene_introns = {}

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        if len(fields) < 9 or fields[2] != "intron":
            continue

        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]

        gene_id = None
        for item in attrs.split(";"):
            item = item.strip()
            if item.startswith("gene_id"):
                gene_id = item.split()[1].strip('"')
                break

        if gene_id is None:
            continue

        gene_introns.setdefault(gene_id, set()).add(
            (chrom, start, end, strand)
        )


rna_supported = set()
protein_supported = set()

for gene_id, introns in gene_introns.items():
    if any(intron in rna_introns for intron in introns):
        rna_supported.add(gene_id)

    if any(intron in protein_introns for intron in introns):
        protein_supported.add(gene_id)


both_supported = rna_supported & protein_supported

n_multiexon = len(gene_introns)

print(f"Multi-exon genes: {n_multiexon}")
print(f"RNA-supported: {len(rna_supported)} ({100 * len(rna_supported) / n_multiexon:.2f}%)")
print(f"Protein-supported: {len(protein_supported)} ({100 * len(protein_supported) / n_multiexon:.2f}%)")
print(f"Supported by both: {len(both_supported)} ({100 * len(both_supported) / n_multiexon:.2f}%)")
