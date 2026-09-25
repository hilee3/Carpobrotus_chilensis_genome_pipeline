import sys


if len(sys.argv) != 3:
    sys.exit(
        "Usage: python 01_codon_fasta_to_axt.py "
        "input.codon.fa output.axt"
    )

infile = sys.argv[1]
outfile = sys.argv[2]

seqs = []
name = None
chunks = []

with open(infile) as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(chunks)))

            name = line[1:].split()[0]
            chunks = []

        else:
            chunks.append(line)

if name is not None:
    seqs.append((name, "".join(chunks)))

if len(seqs) != 2:
    sys.exit(
        f"Expected 2 sequences, found {len(seqs)} in {infile}"
    )

(name1, seq1), (name2, seq2) = seqs

if len(seq1) != len(seq2):
    sys.exit("Aligned sequences have different lengths.")

with open(outfile, "w") as out:
    out.write(f"{name1}\n")
    out.write(f"{seq1}\n")
    out.write(f"{seq2}\n")
