### Divide input fasta file into two output files: one as a reference and one as a query (e.g., for pairwise sequence alignment).

from Bio import SeqIO
import argparse

if __name__ == "__main__":

    # Parse arguments
    argparser = argparse.ArgumentParser(description="Divide input fasta file into a two fasta files: one containing the reference sequence, the other the query sequences.")
    argparser.add_argument("--fasta", help="Input fasta file.")
    argparser.add_argument("--reference_accn", help="Accession number of the reference sequence.")
    argparser.add_argument("--output_ref", help="Output file for the reference sequence.")
    argparser.add_argument("--output_query", help="Output file for the query sequences.")

    args = argparser.parse_args()

    # Read input fasta file
    sequences = list(SeqIO.parse(args.fasta, "fasta"))

    # Divide sequences into reference and query by accession number of the reference sequence
    reference_sequences = []
    query_sequences = []
    for seq in sequences:
        if seq.id in args.reference_accn:
            reference_sequences.append(seq)
        else:
            query_sequences.append(seq)
            
    SeqIO.write(reference_sequences, args.output_ref, "fasta") # Write reference sequences to file
    SeqIO.write(query_sequences, args.output_query, "fasta") # Write query sequences to file