### This script removes duplicate sequences from a multiple sequence alignment.

from Bio import SeqIO
import argparse

if __name__ == "__main__":
        
    argparser = argparse.ArgumentParser(description="Remove duplicate sequences from a multiple sequence alignment.")
    argparser.add_argument("--alignment", help="Multiple sequence alignment in fasta format.")
    argparser.add_argument("--output", help="Output file name.")
        
    args = argparser.parse_args()
        
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    
    # Find sequences that have the same record id
    new_alignment = []
    duplicate_ids = set()
    for record in alignment:
        if record.id in [record.id for record in new_alignment]:
            duplicate_ids.add(record.id)
        else:
            new_alignment.append(record)
    
    # Print duplicate record ids
    if len(duplicate_ids) > 0:
        print("Duplicate sequences found:")
        print("\n".join(duplicate_ids))
    else:
        print("No duplicate sequences found.")
        
    # Write new alignment to a fasta file
    SeqIO.write(new_alignment, args.output, "fasta")
        