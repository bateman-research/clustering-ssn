# Protein sequence similarity network (SSN) and clustering

This repository contains scripts to construct Sequence Similarity Networks (SSNs) of proteins and to cluster them by sequence similarity.

### Software requirements

- NCBI BLAST+: https://blast.ncbi.nlm.nih.gov
- `R` v3.6 with packages `seqinr`, `argparse`, `dplyr`, `tidyr`, `igraph`, `data.table`, `ggplot2` and `tidyr`.

### Scripts

This repository contains three scripts, listed below. 
Each script can be run from the command line with custom input and output options. 
Use the help option (-h/--help) on each script to see all the available options.

```
Rscript script_name.R -h
```

- [blastp_clustering.R](blastp_clustering.R): plot the nucleotide composition and bias profile along the DNA sequence of a gene.
- [pfam_ssn.R](pfam_ssn.R): plot the amino acid composition and bias profile along a protein sequence.
- [pfam_pidmat.R](pfam_pidmat.R): analyse the DNA composition for a subset of genes, including GC content, nucleotide asymmetries and expected vs observed nucleotide frequencies, and generate PCA and dendogram plots.

### Results

Each script takes as input protein sequences in a FASTA format and/or protein alignments in Pfam format, plus other options and parameters, and produces figures and/or results tables.
Results and figures for example input files are provided in the [examples](examples) folder.

### Reference

These scripts were created and used for my PhD research work at the European Bioinformatics Institute (EMBL-EBI). 
More details are provided in my Thesis (insert link when available).

Aleix Lafita - April 2021
