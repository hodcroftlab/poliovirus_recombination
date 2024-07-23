### This script takes a multi-fasta file of coding regions as an input and keeps sequences that start with a gap (as opposed to a start codon)
### Input: a multi-fasta file of coding regions
### Output: a multi-fasta file of coding regions that start with a gap

from Bio import SeqIO

# Read input fasta
fasta = SeqIO.parse(snakemake.input.coding_regions_fasta, "fasta")

def filter_missing_start(fasta):
    filtered_seqs = []
    for record in fasta:
        # Check if the sequence starts with a gap
        if record.seq[0] == "-":
            filtered_seqs.append(record)
    return filtered_seqs

missing_start_fasta = filter_missing_start(fasta)

# Write the filtered sequences to a file
SeqIO.write(missing_start_fasta, snakemake.output.missing_start_fasta, "fasta")