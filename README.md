# Enterovirus C: Phylogenetic & Recombination Analysis

_STILL UNDER CONSTRUCTION_

This repository provides the code for a phylogenetic and recombination analysis of the *Enterovirus C* species, presented in my (unpublished) MSc thesis. 
Whole genome (>6000 bp) enterovirus C sequences (taxid:138950) were downloaded from **[NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Enterovirus%20C,%20taxid:138950)** in May 2024. 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) was used as a workflow management system; the `snakefile` provided in this repository contains the computational pipeline for both the phylogenetic and recombination analysis.
The `scripts` directory contains individual Python and R scripts that are called by the `snakefile`. 

### Phylogenetic Analysis using Nextstrain

Phylogenetic analysis was performed based on the whole genome alignment and alignments of all individual genes as well as the 5' untranslated region using the Nextstrain phylogenetics pipeline.
The interactive trees generated through this code can be interactively explored on the [Nextstrain website](https://nextstrain.org/groups/hodcroftlab/enterovirus/c/VP1).
Refer to the **[Nextstrain publication](https://doi.org/10.1093/bioinformatics/bty407)** and **[Nextstrain documentation](https://docs.nextstrain.org/en/latest/)** for more information on the project. 

### Recombination Analysis

Recombination analysis was performed using a custom similarity plotting approach (inspired by [SimPlot](https://mybiosoftware.com/simplot-3-5-1-sequence-similarity-plotting.html) and [SimPlot++](https://doi.org/10.1093/bioinformatics/btac287)) and the recombination detection method [VirusRecom](https://doi.org/10.1093/bib/bbac513).

... more to come soon!
