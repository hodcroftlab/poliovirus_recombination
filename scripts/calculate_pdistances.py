### This script calculates pairwise distances between all sequences in an alignment fasta file
### The output is a dataframe of pairwise distances (seq1, seq2, pdist_gene1, pdist_gene2)

from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

# p-distance (Hamming distance) is defined as p=nd/n, where n is the number of sites compared and nd is the number of sites that differ between the two sequences

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Calculate pairwise p-distances between all sequences in the given alignment(s).")
    argparser.add_argument("--alignments", nargs="+", help="Input alignment file(s) in fasta format. Multiple alignments can be provided, e.g., to cover multiple genes of the same original sequence.")
    argparser.add_argument("--output", help="Output file (.csv) for a dataframe containing the calculated pairwise distances.")
    
    args = argparser.parse_args()
    
    # Function to calculate pairwise distances between all sequences in an alignment
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
        
    alignment_files = args.alignments
    #alignment_files = ["../data/sequences/alignments/by_cluster/CVA17/CVA17_VP1_aligned.fasta", "../data/sequences/alignments/by_cluster/CVA17/CVA17_3D_aligned.fasta"]
        
    # Make dictionary mapping regions to alignments
    regions_to_alignment = {}
    for alignment_file in alignment_files:
        region = alignment_file.split("/")[-1].split("_")[-2]
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        regions_to_alignment[region] = alignment
    
    # Calculate pairwise distances for each alignment
    results_dict = {}
    for region, alignment in regions_to_alignment.items():
        print(f"\nCalculating p-distances for {region} ...")
        pdist_df = calculate_pairwise_distances(alignment)
        results_dict[region] = pdist_df
        
        
    # Merge the results from all alignments into a single dataframe
    print("\nMerging results ...")
    merged_df = None
    for region, pdist_df in results_dict.items():
        # Rename the pdist column to the region name
        pdist_df = pdist_df.rename(columns={"pdist": f"{region}"})
        if merged_df is None:
            merged_df = pdist_df
        else:
            merged_df = merged_df.merge(pdist_df, on=["seq1", "seq2"], how="outer")
            
    # Write the merged dataframe to a CSV file
    merged_df.to_csv(args.output, index=False)

    
    