### This script builds a consensus sequence from the aligned sequences in the input file.

from Bio import SeqIO
import argparse
import numpy as np
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build consensus sequence from the aligned sequences in the input file.")
    parser.add_argument("--input", nargs = "+", help="Input alignment file(s) (.fasta).")
    parser.add_argument("--output", help="Output file to write the consensus sequence (.fasta).")
    args = parser.parse_args()
    
    # Read the aligned sequences
    alignment_dict = {}
    for aln in args.input:
        name = aln.split("/")[-1].split(".")[0]
        alignment_dict[name] = list(SeqIO.parse(aln, "fasta"))

    consensus_dict = {}
    for name, alignment in alignment_dict.items():
        
        # Get number of sequences
        num_sequences = len(alignment)
        
        # Convert alignment to a numpy array
        alignment_array = np.array([list(str(record.seq)) for record in alignment])
        
        # Find the most common base at each position to build the consensus sequence
        consensus_sequence = ""
        for i in range(alignment_array.shape[1]):
            bases, counts = np.unique(alignment_array[:, i], return_counts=True)
            consensus_sequence += bases[np.argmax(counts)]
            
        # Store the consensus sequence
        consensus_dict[name] = (consensus_sequence, num_sequences)
        
    # Write the consensus sequences to the output file (as fasta)
    with open(args.output, "w") as output_handle:
        for name, (consensus_sequence, num_sequences) in consensus_dict.items():
            name = name.replace("all", "EVC").replace("non_polio", "non-polio").replace("_", "/")
            output_handle.write(f">{name} consensus sequence, built from {num_sequences} {name} sequences\n")
            output_handle.write(f"{consensus_sequence}\n")
    
