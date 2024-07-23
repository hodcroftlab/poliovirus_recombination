### This script removes trailing nucleotides (not part of a codon) from the start and end of a nucleotide alignment (note: this removes all gaps and unaligns the sequences).

from Bio import SeqIO
import argparse
import json

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Remove trailing nucleotides (that are not part of a codon) from the start and end of a nucleotide alignment (note: this removes all gaps and unaligns the sequences).")
    
    argparser.add_argument("--fasta", help="Input fasta file")
    argparser.add_argument("--trailing_starts", help="Dictionary (json) of sequence IDs and number of nucleotides that need to be removed from the start of these sequences")
    argparser.add_argument("--output", help="Output fasta file")

    args = argparser.parse_args()

    # Read input fasta
    fasta = SeqIO.parse(args.fasta, "fasta")

    # Read dictionary
    with open(args.trailing_starts) as f:
        trailing_starts = json.load(f)

    # Remove trailing nucleotides
    cleaned_seqs = []
    for record in fasta:
        # Remove gaps
        record.seq = record.seq.replace("-", "")
        if record.id in trailing_starts: # Remove nucleotides from the start    
            record.seq = record.seq[trailing_starts[record.id]:]
        record.seq = record.seq[:len(record.seq) - len(record.seq) % 3] # Remove nucleotides from the end by checking if the sequence is a multiple of 3
        cleaned_seqs.append(record)
        
    # Write the cleaned sequences to a file
    SeqIO.write(cleaned_seqs, args.output, "fasta")
