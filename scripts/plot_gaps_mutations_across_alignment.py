### This script counts the number of gaps per amino acid and mutations compared to the reference sequence in a multiple sequence nucleotide alignment in fasta format.

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Count the number of gaps per amino acid / codon position in a multiple sequence nucleotide alignment in fasta format. Assumes in-frame alignment (first reading frame)")
    parser.add_argument("--alignments", help="Multiple sequence alignment(s) in fasta format.", nargs="+")
    parser.add_argument("--reference_accn", help="Accession number of the reference sequence in the alignment.")
    parser.add_argument("--gene_annotation", help="Gene annotation file in .csv format (expected columns: feature, start, end).")
    parser.add_argument("--output_csv", help="Output directory for .csv files containing the number and proportion of gaps and mutations per amino acid position for each alignment.")
    parser.add_argument("--output_plot", help="Output directory for .png files containing a plot of the distribution of gaps and mutations per amino acid position for each alignment.")
    args = parser.parse_args()

    # Iterate over the alignment files
    for alignment_file in args.alignments:
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        sequence_length = len(alignment[0].seq) 
        sequence_number = len(alignment)
        
        # Extract name of the alignment file
        alignment_name = alignment_file.split("/")[-1]
        # Remove .fasta extension
        if alignment_name.endswith(".fasta"):
            alignment_name = alignment_name[:-6]
        
        print(f"Alignment {alignment_name} contains {sequence_number} sequences of length {sequence_length} nucleotides.")
        
        # # Use BioPython to translate the sequences to amino acids
        # translated_seqs = []
        
        # for record in alignment:
        #     # Convert incomplete codons to gaps
        #     sequence = str(record.seq)
        #     new_sequence = []
        #     for i in range(0, len(sequence), 3):
        #         codon = sequence[i:i+3]
        #         if "-" in codon and codon.count("-") < 3:
        #             new_sequence.append("---")
        #         else:
        #             new_sequence.append(codon)
        #     record.seq = Seq("".join(new_sequence))
                    
        #     # Translate the sequence
        #     record.seq = record.seq.translate()
        #     translated_seqs.append(record)
            
        # Get the reference sequence
        reference_seq = None
        for record in alignment:
            if record.id == args.reference_accn:
                reference_seq = record.seq
                break

        if reference_seq is None:
            raise ValueError(f"Reference sequence with ID {args.reference_accn} not found in alignment.")

        # Initialize a dictionary to store the number of gaps and mutations per amino acid position
        position_dict = {}
        
        # Iterate over the sequences in the alignment
        for record in alignment:
            sequence = str(record.seq)
            for j in range(len(sequence)):
                nt = sequence[j]
                if j not in position_dict:
                    position_dict[j] = {"gaps": 0, "mutations": 0}
                if nt == "-":
                    position_dict[j]["gaps"] += 1
                elif nt != reference_seq[j]:
                    position_dict[j]["mutations"] += 1
                    
        # Convert the dictionary to a dataframe with three columns: aa_position, gaps, mutations
        position_df = pd.DataFrame(position_dict).T
        # Move the index to a column
        position_df.reset_index(inplace=True)
        
        # Convert the number of gaps and mutations to proportions
        position_df["gaps_proportion"] = round(position_df["gaps"] / sequence_number, 4)
        position_df["mutations_proportion"] = round(position_df["mutations"] / sequence_number, 4)
        
        # Save dataframe to a .csv file
        output_csv = args.output_csv + "/" + alignment_name + ".csv"
        position_df.to_csv(output_csv, index=False)
        
        
        ## Plot distribution of gaps and mutations per amino acid position in the alignment
        # Density plot, with gaps and mutations as separate lines using seaborn
        # Add gene annotation to the plot (as bars)
        
        # Set figure style and dimensions
        sns.set_style("whitegrid")
        plt.figure(figsize=(24, 6))
        
        # Make lineplots for gaps and mutations
        sns.lineplot(x="index", y="mutations_proportion", data=position_df, color="grey", label="Mutations", alpha=0.5, linewidth=0.5)
        sns.lineplot(x="index", y="gaps_proportion", data=position_df, color="black", label="Gaps", linewidth=1.2)
        
        # Adjust axis limits and labels
        plt.ylim(-0.1, 1) # Adjust to leave space below the plot for gene annotations
        plt.ylabel("Proportion of sequences")
        plt.xlim(0, sequence_length-1)
        plt.xlabel("Nucleotide position")
        # Adjust xticks to show move in steps of 500
        plt.xticks(range(0, sequence_length, 500))
        
        # Read the gene annotation file
        gene_annotation = pd.read_csv(args.gene_annotation, sep="\t")
        # Subtract 1 from all coordinates to convert from 1-based to 0-based
        gene_annotation["start"] -= 1
        gene_annotation["end"] -= 1
        
        # Get colors for the gene annotations
        colors_annotation = sns.color_palette("tab20", n_colors=len(gene_annotation))
        
        # Add gene annotations as bars below the plot
        for i in range(len(gene_annotation)):
            start = gene_annotation.iloc[i]["start"]
            end = gene_annotation.iloc[i]["end"]
            gene_name = gene_annotation.iloc[i]["feature"]
            rect = patches.Rectangle((start, -0.1), end-start, 0.1, linewidth=1, edgecolor=colors_annotation[i], facecolor=colors_annotation[i], alpha=0.5, label=gene_name)
            plt.gca().add_patch(rect)
            plt.text((start + end) / 2, -0.05, gene_name, ha="center", va="center", fontsize=8, color="black", rotation=0)
        
        # Save plot to a .png file
        plt.tight_layout() # Tight layout
        output_plot = args.output_plot + "/" + alignment_name + ".png"
        plt.savefig(output_plot)
        
    
    # # Lastly, make a combined plot of all alignments in the same figure for the gaps only
    # plt.figure(figsize=(24, 6))
    # # Use sns colorblind friendly palette
    # sns.set_palette("colorblind")
    # # Make lineplots for gaps
    # for alignment_file in args.alignments:
    #     alignment_name = alignment_file.split("/")[-1]
    #     # Remove .fasta extension
    #     if alignment_name.endswith(".fasta"):
    #         alignment_name = alignment_name[:-6]
    #     position_df = pd.read_csv(args.output_csv + "/" + alignment_name + ".csv")
    #     sns.lineplot(x="index", y="gaps_proportion", data=position_df, label=alignment_name, linewidth=1.5, alpha=0.8)
        
    # # Adjust axis limits and labels
    # plt.ylim(-0.1, 1) # Adjust to leave space below the plot for gene annotations
    # plt.ylabel("Proportion of sequences")
    # plt.xlim(0, aa_sequence_length-1)
    # plt.xlabel("Amino acid position")
    # # Adjust xticks to show move in steps of 200
    # plt.xticks(range(0, aa_sequence_length, 200))
    # plt.tight_layout()
    
    # # Again, add gene annotations as bars below the plot
    # for i in range(len(gene_annotation)):
    #     start = gene_annotation.iloc[i]["start"]
    #     end = gene_annotation.iloc[i]["end"]
    #     gene_name = gene_annotation.iloc[i]["feature"]
    #     rect = patches.Rectangle((start, -0.1), end-start, 0.1, linewidth=1, edgecolor=colors_annotation[i], facecolor=colors_annotation[i], alpha=0.5, label=gene_name)
    #     plt.gca().add_patch(rect)
    #     plt.text((start + end) / 2, -0.05, gene_name, ha="center", va="center", fontsize=8, color="black", rotation=0)
        
    # plt.savefig(args.output_plot + "/gaps_proportion_combined.png")
    

        