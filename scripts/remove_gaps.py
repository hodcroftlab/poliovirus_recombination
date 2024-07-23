### This script removes gap characters ("-") from all sequences in a fasta file

from Bio import SeqIO
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="remove_gaps", description="Remove gap characters ('-') from all sequences in a fasta file.")
    parser.add_argument("--input", help="Path to the input fasta file.", type=str)
    parser.add_argument("--output", help="Path to the output fasta file.", type=str)
    args = parser.parse_args()
    
    input_fasta = args.input
    output_fasta = args.output

    with open(input_fasta, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    for record in records:
        record.seq = record.seq.replace("-", "")
        
    # Write the records to a new fasta file
    with open(output_fasta, "w") as handle:
        SeqIO.write(records, handle, "fasta")