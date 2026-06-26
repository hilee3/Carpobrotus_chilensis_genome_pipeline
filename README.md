# C. chilensis Genome Assembly and Annotation Pipeline

This repository contains the scripts used to assemble, evaluate, and annotate the *Carpobrotus chilensis* reference genome from PacBio HiFi sequencing and RNA-seq data.

## Directory Structure

* **01_Assembly**

  * Initial genome assembly with Hifiasm and haplotig purging using purge_dups.

* **02_Contamination_filtering**

  * Identification and removal of organelle-derived and other contaminant contigs.

* **03_QC**

  * Assembly quality assessment using BUSCO, QUAST, Meryl, and Merqury.

* **04_Repeat_annotation**

  * Repeat annotation and genome soft-masking using EDTA and RepeatMasker.

* **05_Gene_annotation**

  * RNA-seq alignment, gene prediction with BRAKER3, and annotation processing.

* **06_Functional_annotation**

  * Functional annotation using eggNOG-mapper.

## Requirements

Software requirements, input files, and computational resource specifications are listed in the header of each script.

## Repository Organization

The scripts are organized according to the order in which they were used during genome assembly and annotation. Each directory corresponds to a major stage of the analysis.

## Notes

## Notes

The scripts are written as SLURM batch scripts for execution on an HPC system. Input and output paths, computational resources, and scheduler directives may require modification for other computing systems.

## Citation

If you use this repository, please cite the associated manuscript.

## License

This repository is released under the MIT License.