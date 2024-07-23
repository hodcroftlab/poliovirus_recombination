### Script to keep only the first n sequences of a given fasta file
from Bio import SeqIO

# Define function to filter fasta sequences by length using BioPython
def get_first_n_seqs(fasta, n_seqs):
    sequences_to_keep = []
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences_to_keep.append(record)
            if len(sequences_to_keep) >= n_seqs:
                break
    return sequences_to_keep

# Apply function on snakemake with snakemake inputs/parameters
sequences_to_keep = get_first_n_seqs(fasta=snakemake.input.fasta, n_seqs=snakemake.params.n_seqs)

# Write filtered sequences to output files            
SeqIO.write(sequences_to_keep, snakemake.output.sliced_fasta, "fasta")