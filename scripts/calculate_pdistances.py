### This script calculates pairwise distances between all sequences in an alignment fasta file
### The output is a dataframe of pairwise distances (seq1, seq2, pdistance)

from Bio import SeqIO
import pandas as pd
import numpy as np


# p-distance is defined as p=nd/n, where n is the number of sites compared and nd is the number of sites that differ between the two sequences

def calculate_pdistances(alignment_path):
    # Read alignment file into a list of SeqRecord objects
    alignment = list(SeqIO.parse(alignment_path, "fasta"))
    
    # Filter out sequences with excessive ambiguous characters
    valid_chars = set("ACGTacgt-")
    threshold = 0.1 # Maximum fraction of non-valid characters allowed
    alignment = [rec for rec in alignment if sum(c not in valid_chars for c in rec.seq) / len(rec.seq) < threshold]
    
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
            
            # Calculate the number of differing positions (including gaps) --> this also includes N=N as a match, but that should be very rare
            nd = np.sum(seq1 != seq2)
            
            # Calculate the p-distance
            p = nd / seq_len
            
            # Append the result as a tuple to the results list
            results.append((alignment[i].id, alignment[j].id, p))
            
    # Convert the list of results to a dataframe
    results_df = pd.DataFrame(results, columns=["seq1", "seq2", "pdistance"])
    
    return results_df
            

# Run the function on the specified path to the alignment
pdistances_vp1 = calculate_pdistances(snakemake.input.alignment_vp1).rename(columns={"pdistance": "dist_vp1"})
pdistances_3dpol = calculate_pdistances(snakemake.input.alignment_3dpol).rename(columns={"pdistance": "dist_3dpol"})
pdistances_5utr = calculate_pdistances(snakemake.input.alignment_5utr).rename(columns={"pdistance": "dist_5utr"})

# Merge the three dataframes
pdistances = pdistances_vp1.merge(pdistances_3dpol, on=["seq1", "seq2"], how="outer").merge(pdistances_5utr, on=["seq1", "seq2"], how="outer")

# Sort by seq1 and seq2
pdistances = pdistances.sort_values(["seq1", "seq2"])

# Write the merged dataframe to the output file
pdistances.to_csv(snakemake.output.pdistances_csv, index=False)