### This script removes insertions in a multiple sequece alignment relative to a reference sequence.

from Bio import SeqIO
from Bio.Seq import Seq
import argparse

if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description="Remove insertions relative to a reference sequence in a multiple sequence alignment.")
    argparser.add_argument("--alignment", help="Multiple sequence alignment in fasta format.")
    argparser.add_argument("--reference", help="Identifier for the reference sequence in the multiple sequence alignment (sequence relative to which gaps will be removed).")
    argparser.add_argument("--output", help="Output file name.")

    args = argparser.parse_args()

    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    reference = args.reference

    # Get reference sequence
    reference_sequence = None
    for record in alignment:
        if record.id == reference:
            reference_sequence = record.seq
            break

    if reference_sequence is None:
        raise ValueError(f"Reference sequence with ID {reference} not found in the alignment.")
            
    # Remove insertions relative to the reference sequence and write them to a new fasta file
    new_alignment = []
    for record in alignment:
        new_sequence = Seq("")
        for i in range(len(reference_sequence)):
            if reference_sequence[i] != "-":
                new_sequence += record.seq[i]
        record.seq = new_sequence
        new_alignment.append(record)
            
    SeqIO.write(new_alignment, args.output, "fasta")


