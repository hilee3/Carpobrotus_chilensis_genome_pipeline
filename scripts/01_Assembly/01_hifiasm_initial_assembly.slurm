#!/bin/bash

#SBATCH --job-name=Hifiasm_chilensis
#SBATCH --time=10-01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=Hifiasm_chilensis.out
#SBATCH --error=Hifiasm_chilensis.err

# De novo assembly of the C. chilensis PacBio HiFi reads using hifiasm.
# Replace the input path below with the local path to the HiFi FASTQ file.

READS="/path/to/m84066_241205_191147_s1.fastq"
PREFIX="chilensis_assembly_l0_031025.asm"

hifiasm \
    -o "${PREFIX}" \
    -t 32 \
    "${READS}"