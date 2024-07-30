### This script calculates pairwise distances between all sequences in an alignment fasta file
### The output is a dataframe of pairwise distances (seq1, seq2, pdist_gene1, pdist_gene2)

from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

# p-distance is defined as p=nd/n, where n is the number of sites compared and nd is the number of sites that differ between the two sequences

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Calculate pairwise p-distances between all sequences in the given alignment(s).")
    argparser.add_argument("--alignments", help="Input alignment file(s) in fasta format. Multiple alignments can be provided to cover multiple genes of the same original sequence.", nargs="+")
    argparser.add_argument("--output", help="Output file (.csv) for a dataframe containing the calculated pairwise distances.")
    
    args = argparser.parse_args()
    
    # Read alignment files into a dictionary
    alignments = {}
    for alignment_file in args.alignments:
        name = alignment_file.split("/")[-1].split("_")[0]
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        alignments[name] = alignment
    
    def calculate_pairwise_distances(alignment):
        n_seqs = len(alignment)
        print("Number of sequences: ", n_seqs)
        seq_len = len(alignment[0].seq)
        print("Sequence length: ", seq_len)
        
        # Convert the sequences to a NumPy array for faster element-wise calculations
        alignment_array = np.array([list(rec.seq) for rec in alignment])
        
        # Initialize empty list to collect results
        results = []
        
        for i in range(n_seqs):
            seq1 = alignment_array[i]
            for j in range(i, n_seqs): # Start j from i to avoid duplicates
                
                seq2 = alignment_array[j]
                
                # Mask all positions where either of the two sequences has a gap or an N
                valid_positions = (seq1 != "-") & (seq1 != "N") & (seq2 != "-") & (seq2 != "N")
                seq1_valid = seq1[valid_positions]
                seq2_valid = seq2[valid_positions]
                
                # Update the sequence length
                seq_len_valid = len(seq1_valid)
                
                # Now calculate the number of differing positions
                nd = np.sum(seq1_valid != seq2_valid)
                
                # Calculate the p-distance
                pdist = nd / seq_len_valid
                pdist = round(pdist, 4) # Round to 4 decimal places
                
                # Append the result as a tuple to the results list
                results.append((alignment[i].id, alignment[j].id, pdist))
                
        # Convert the list of results to a dataframe
        results_df = pd.DataFrame(results, columns=["seq1", "seq2", "pdist"])
        
        return results_df
        
    # Iterate over alignments in the dictionary and calculate pairwise distances
    results_dict = {}
    for name, alignment in alignments.items():
        print(f"Calculating pairwise distances for {name}...")
        pdist_df = calculate_pairwise_distances(alignment)
        results_dict[name] = pdist_df
    
    # Merge the results from all alignments by seq1 and seq2 columns; add a suffix to the p-distance columns to distinguish between different genes
    merged_df = None
    for name, df in results_dict.items():
        # Rename the pdistance column to include the gene name
        df = df.rename(columns={"pdist": f"pdist_{name}"})
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=["seq1", "seq2"], how="outer", )
        
        
    
    # Write the merged dataframe to a CSV file
    merged_df.to_csv(args.output, index=False)