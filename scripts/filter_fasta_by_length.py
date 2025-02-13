### Script to filter fasta sequences by length using BioPython
# Input: fasta file with multiple sequences; minimum length of sequences to keep
# Output: filtered fasta file with sequences longer than the minimum length

from Bio import SeqIO
import argparse
import pandas as pd

# Define argparse arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("--fasta", help="Input fasta file.")
argparser.add_argument("--metadata", help="Input metadata file.")
argparser.add_argument("--min_length", type=int, help="Minimum length of sequences to keep.")
argparser.add_argument("--output_fasta", help="Output filtered fasta file.")
argparser.add_argument("--output_metadata", help="Output filtered metadata file.")

args = argparser.parse_args()

# Filter fasta sequences by length using BioPython
sequences_to_keep = []
with open(args.fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) >= args.min_length:
            sequences_to_keep.append(record)

print("Number of sequences after filtering: ", len(sequences_to_keep))

# Filter metadata based on sequences that were kept
metadata = pd.read_csv(args.metadata) # Read metadata file
metadata_filtered = metadata[metadata["Accession"].isin([seq.id for seq in sequences_to_keep])] # Filter metadata based on sequences that were kept

# Write filtered sequences and metadata to output files            
SeqIO.write(sequences_to_keep, args.output_fasta, "fasta")
metadata_filtered.to_csv(args.output_metadata, index=False)