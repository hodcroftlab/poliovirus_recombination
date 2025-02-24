# Enterovirus C: Phylogenetic & Recombination Analysis

_STILL UNDER CONSTRUCTION_

This repository provides the code for a phylogenetic and recombination analysis of the *Enterovirus C* species, presented in my (unpublished) MSc thesis. 

Whole genome (>6000 bp) enterovirus C sequences (taxid:138950) were downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Enterovirus%20C,%20taxid:138950) in May 2024; a total of 1940 sequences EV-C sequences were retained for the analysis.

[Snakemake](https://snakemake.readthedocs.io/en/stable/) was used as a workflow management system; the `Snakefile` provided in this repository contains the computational pipeline for both the phylogenetic and recombination analysis.
The `scripts` directory contains individual Python and R scripts that are called by the `Snakefile`. 


### Phylogenetic Analysis using Nextstrain

Phylogenetic analysis was performed based on the whole genome alignment and alignments of all individual genes as well as the 5' untranslated region using the **Nextstrain** phylogenetics pipeline. The interactive trees generated through this code can be interactively explored on the [Nextstrain website](https://nextstrain.org/groups/hodcroftlab/enterovirus/c/VP1).

To install the Nextstrain environment, follow [these instructions](https://docs.nextstrain.org/en/latest/install.html). Once the Nextstrain environment has been set up and activated, phylogenetic analysis can be performed by executing `snakemake --cores 1 export_all`. Alternatively, the rules specified in the `Snakefile` can be executed individually in a step-wise manner. The generated trees (in JSON format) can be visualized using the `auspice view` command (not included in the `Snakefile`; use `auspice view -h` for help). 

Refer to the [Nextstrain publication](https://doi.org/10.1093/bioinformatics/bty407) and [Nextstrain documentation](https://docs.nextstrain.org/en/latest/) for more information on the project. 


### Recombination Analysis

Recombination analysis was performed using a custom similarity plotting approach (inspired by [SimPlot](https://mybiosoftware.com/simplot-3-5-1-sequence-similarity-plotting.html) and [SimPlot++](https://doi.org/10.1093/bioinformatics/btac287)) and the recombination detection method [VirusRecom](https://doi.org/10.1093/bib/bbac513). The code to generate the similarity plots is provided in the `scripts/custom_simplots.py` and `scripts/custom_simplots_extended.py` scripts.

The similarity plots and VirusRecom results for all 1940 sequences can be downloaded [here](https://drive.google.com/drive/folders/1UzJHC-0qk4CdhvUEeS8Zg--sN9TxLO65?usp=sharing) for exploration.

... more to come soon!
