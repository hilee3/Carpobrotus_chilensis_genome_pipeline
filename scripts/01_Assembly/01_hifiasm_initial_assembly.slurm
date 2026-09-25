#!/bin/bash

#SBATCH --job-name=Hifiasm_chilensis
#SBATCH --time=10-01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=Hifiasm_chilensis.out
#SBATCH --error=Hifiasm_chilensis.err

# De novo assembly of C. chilensis PacBio HiFi reads using hifiasm.
# Replace the input path below with the local path to the HiFi FASTQ file.

READS="/path/to/pacbio_hifi_reads.fastq"
PREFIX="chilensis_assembly.asm"

hifiasm \
    -o "${PREFIX}" \
    -t 32 \
    "${READS}"

# Convert the primary contig assembly from GFA to FASTA.
awk '/^S/{print ">"$2; print $3}' \
    "${PREFIX}.bp.p_ctg.gfa" \
    > "${PREFIX}.bp.p_ctg.fa"