### Script to filter fasta sequences by length using BioPython
# Input: fasta file with multiple sequences; minimum length of sequences to keep
# Output: filtered fasta file with sequences longer than the minimum length

from Bio import SeqIO

# Define function to filter fasta sequences by length using BioPython
def filter_fasta_by_length(fasta, min_length):
    sequences_to_keep = []
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record.seq) >= min_length:
                sequences_to_keep.append(record)
    return sequences_to_keep

# Apply function on snakemake with snakemake inputs/parameters
sequences_to_keep = filter_fasta_by_length(fasta=snakemake.input.fasta, min_length=snakemake.params.min_length)
print("Number of sequences after filtering: ", len(sequences_to_keep))

# Write filtered sequences to output files            
SeqIO.write(sequences_to_keep, snakemake.output.filtered_fasta, "fasta")