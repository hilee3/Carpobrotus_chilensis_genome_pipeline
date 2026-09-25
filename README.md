# *Carpobrotus chilensis* Genome Assembly, Annotation, and Comparative Analysis

This repository contains the scripts used for assembly, quality assessment, and annotation of the *Carpobrotus chilensis* draft genome using PacBio HiFi sequencing and RNA-seq data. It also includes scripts used for reannotation of the published *C. edulis* genome, orthology inference, pairwise dN/dS analysis, and GO enrichment analysis.

## Repository structure

* **01_Assembly**
  * Initial genome assembly with Hifiasm, haplotig purging with purge_dups, and chloroplast genome assembly with OATK.

* **02_Contamination_filtering**
  * Identification and removal of organelle-derived contigs, additional chloroplast screening, contaminant screening with Tiara, and correction of adaptor-like sequences identified by FCS-GX.

* **03_Assembly_QC**
  * Assembly quality assessment using BUSCO, QUAST, and Merqury; HiFi read mapping and alignment-quality assessment; sequencing-depth and coverage-breadth calculations; and BlobToolKit snailplot generation.

* **04_Repeat_annotation**
  * Repeat annotation with EDTA, genome soft-masking with RepeatMasker, and assessment of repeat-space assembly quality using the LTR Assembly Index (LAI).

* **05_Gene_annotation**
  * RNA-seq alignment with HISAT2, gene prediction with BRAKER3, longest-isoform selection with AGAT, protein and CDS extraction with gffread, and annotation statistics.
  * Additional annotation-quality assessments include protein-mode BUSCO analysis, quantification of RNA-seq and protein evidence support for multi-exon gene models, overlap between predicted gene models and annotated transposable elements, and screening of predicted proteins for TE-associated Pfam annotations.

* **06_Functional_annotation**
  * Functional annotation of predicted *C. chilensis* proteins using eggNOG-mapper.

* **07_C_edulis_reannotation**
  * Reannotation of the published *C. edulis* nuclear genome using EDTA and RepeatMasker for repeat annotation and soft-masking, two publicly available RNA-seq datasets aligned with HISAT2, and BRAKER3 gene prediction using Viridiplantae protein evidence from OrthoDB v12.
  * Representative longest isoforms were selected with AGAT, protein and CDS sequences were extracted with gffread, and the resulting annotation was evaluated using annotation statistics, protein-mode BUSCO, and eggNOG-mapper.

* **08_Orthology**
  * Orthology inference with OrthoFinder and extraction of one-to-one single-copy ortholog pairs between *C. chilensis* and *C. edulis*.

* **09_dNdS**
  * Protein and codon alignment processing, pairwise dN/dS estimation with KaKs_Calculator, quality filtering, and visualization of dS and dN/dS distributions.

* **10_GO_enrichment**
  * Preparation of gene-to-GO mappings and study/background gene sets, GO enrichment analysis with topGO, multiple-testing correction, and visualization of enrichment results.

## Notes

Scripts are provided in the form used for the analyses, with local file paths replaced by generic placeholders where appropriate. Most computationally intensive analyses were run as SLURM batch jobs on an HPC system. Additional Bash, Python, and R scripts were used for data processing, filtering, statistical analysis, and visualization.

Input and output paths, computational resources, software environments, and scheduler directives may require modification for use on other computing systems.

## Data availability

Data generated in this study are available under NCBI BioProject **PRJNA1468704**.

- Genome assembly: **JBZCIG000000000**
- PacBio HiFi reads: **SRR38808173**
- RNA-seq reads: **SRR38839563**

Additional data and analysis resources are available through Zenodo:

**10.5281/zenodo.20687645**

These resources include genome annotation files, repeat annotations and the repeat library, chloroplast genome resources, *Carpobrotus edulis* reannotation files, orthogroup assignments, codon alignments, dN/dS results, and GO enrichment input and output tables.

The *C. edulis* genome assembly used for comparative analyses is available under accession **GCA_965788445.1**. Public RNA-seq datasets used for *C. edulis* gene prediction are available from ENA under accessions **ERR15315142** and **ERR15315144**.

## Citation

If you use these scripts, please cite the associated manuscript.

## License

This repository is released under the MIT License.