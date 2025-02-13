### This script takes alignments for the UTRs and coding sequences and assembles them into full-length genomes by concatenating the sequences according to their names.

from Bio import SeqIO
import argparse

argparser = argparse.ArgumentParser(description="Assemble full-length genomes from UTR and coding sequence alignments.")
argparser.add_argument("--coding_region", help="Input alignment file for the coding sequences in fasta format.")
argparser.add_argument("--five_prime", help="Input alignment file for the 5' UTR sequences in fasta format.")
argparser.add_argument("--three_prime", help="Input alignment file for the 3' UTR sequences in fasta format.")
argparser.add_argument("--output", help="Output file for the full-length genomes in fasta format.")

args = argparser.parse_args()

# Read the alignments
coding_alignment = list(SeqIO.parse(args.coding_region, "fasta"))
five_prime_alignment = list(SeqIO.parse(args.five_prime, "fasta"))
three_prime_alignment = list(SeqIO.parse(args.three_prime, "fasta"))

# Assemble full-length genomes
full_genomes = []
for coding_record in coding_alignment:
    coding_name = coding_record.id
    five_prime_record = [record for record in five_prime_alignment if record.id == coding_name][0]
    three_prime_record = [record for record in three_prime_alignment if record.id == coding_name][0]
    
    full_genome = coding_record
    full_genome.seq = five_prime_record.seq + coding_record.seq + three_prime_record.seq
    full_genomes.append(full_genome)
    
# Write the full-length genomes to a file
SeqIO.write(full_genomes, args.output, "fasta")