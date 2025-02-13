### This script takes the masked alignment and gene annotation as input and outputs genewise and whole-genome masking statistics (percentage of masked positions)
# Masked positions are annotated as "+" in the alignment

from Bio import SeqIO
import pandas as pd
import argparse
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Calculate percentage of masked positions in a multiple sequence alignment after applying masking algorithm.")
    argparser.add_argument("--alignment", help="Input alignment file to be evaluated.")
    argparser.add_argument("--mask_char", help="Character used for masked positions.")
    argparser.add_argument("--gene_annotation", help="Input gene annotation file (CSV or TSV) containing start and end positions of genes.")
    argparser.add_argument("--outcsv", help="Output file (CSV) containing the percentage of masked positions in each gene.")
    argparser.add_argument("--outplot", help="Output file (PNG) containing a bar plot of the percentage of masked positions in each gene.")
    args = argparser.parse_args()
    
    # Read alignment
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    #alignment = list(SeqIO.parse("../data/sequences/alignments/coding_region/nextalign/nextalign_coding_region_masked_ws15.fasta", "fasta"))
    
    # Since the same positions are masked in all sequences, we can calculate the percentage of masked positions in one sequence
    sequence = alignment[0].seq
    
    # Read gene annotation
    if args.gene_annotation.endswith(".csv"):
        gene_annotations = pd.read_csv(args.gene_annotation)
    elif args.gene_annotation.endswith(".tsv"):
        gene_annotations = pd.read_csv(args.gene_annotation, sep="\t")
    else:  
        print("Gene annotation file must be in CSV or TSV format.")
        
    #gene_annotations = pd.read_csv("../data/gene_annotation/annotation_coding_region.csv", sep="\t")
    
    # Split sequence into sub-sequences for each gene
    gene_sequences = {}
    for feature in gene_annotations["feature"]:
        start = gene_annotations[gene_annotations["feature"] == feature]["start"].values[0]-1
        end = gene_annotations[gene_annotations["feature"] == feature]["end"].values[0]
        sub_sequence = sequence[start:end]
        gene_sequences[feature] = sub_sequence
        
    # Add entire polyprotein to gene_sequences
    gene_sequences["Total"] = sequence
    
    # Calculate percentage of masked positions in each gene and for the whole genome
    mask_char = args.mask_char
    #mask_char = "+"
    gene_masking_stats = {}
    for gene in gene_sequences:
        masked_positions = gene_sequences[gene].count(mask_char)
        total_positions = len(gene_sequences[gene])
        percentage_masked = (masked_positions / total_positions) * 100
        gene_masking_stats[gene] = percentage_masked
    
    # Save results to CSV
    gene_masking_stats_df = pd.DataFrame.from_dict(gene_masking_stats, orient="index", columns=["percentage_masked"])
    gene_masking_stats_df.to_csv(args.outcsv)
    
    ### Plot results
    # Set font
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 18
    
    # Assign colors to genes
    tab20_colors = plt.colormaps["tab20"].colors
    n_genes = len(gene_sequences)  # Number of genes
    gene_colors = [tab20_colors[i % 20] for i in range(n_genes)]  # Cycle through the 20 colors
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot bar chart with custom colors and reduced bar gaps
    bars = ax.bar(gene_masking_stats_df.index, gene_masking_stats_df["percentage_masked"], color=gene_colors, width=0.8) 
    
    # Add labels on top of each bar
    ax.set_ylim(0, 1.15 * gene_masking_stats_df["percentage_masked"].max()) # Extend the y-axis limits to make space for the labels
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,  # Center the text
            height,  # y-coordinate of the text
            f"{height:.1f}",  # Format the value
            ha="center",
            va="bottom"
        )
    # Rotate x-axis labels
    plt.xticks(rotation=90)
    
    ax.set_ylabel("Masked positions (%)")    
    plt.tight_layout()
    plt.savefig(args.outplot)