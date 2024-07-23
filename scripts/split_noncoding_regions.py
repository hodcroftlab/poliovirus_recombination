
### This script takes a fasta file of a full genome sequence alignment as an input and cuts-out the 5'UTR and 3'UTR
### Input: a fasta file of the full genome sequence alignment (KEEPING insertions)
### Outputs: 1) a fasta file of the 5'UTR; 2) for the 3'UTR sequences; 3) for the coding regions

from Bio import SeqIO
import argparse

if __name__ == "__main__":
    
    # Parse arguments
    argparser = argparse.ArgumentParser(description="Split the sequences into 5'UTR, 3'UTR, and coding regions.")
    argparser.add_argument("--alignment", help="Aligned fasta file of full genome sequences.")
    argparser.add_argument("--cutoff_5UTR", help="Position of the first base of the coding region.")
    argparser.add_argument("--cutoff_3UTR", help="Position of the last base of the coding region.")
    argparser.add_argument("--output_coding_region", help="Output file for the coding region.")
    argparser.add_argument("--output_5UTR", help="Output file for the 5'UTR.")
    argparser.add_argument("--output_3UTR", help="Output file for the 3'UTR.")
    
    args = argparser.parse_args()
    
    # Read input fasta
    fasta = SeqIO.parse(args.alignment, "fasta")

    # Define the 5'UTR and 3'UTR regions
    cutoff_5utr = int(args.cutoff_5UTR) # first base of the coding region
    cutoff_3utr = int(args.cutoff_3UTR) # last base of the coding region

    # Split the sequences
    five_utr_seqs = []
    three_utr_seqs = []
    coding_region_seqs = []
        
    for record in fasta:
        # 5'UTR
        five_utr = record[:cutoff_5utr-1]
        five_utr_seqs.append(five_utr)
            
        # 3'UTR
        three_utr = record[cutoff_3utr:]
        three_utr_seqs.append(three_utr)
            
        # Coding region
        coding_region = record[cutoff_5utr-1:cutoff_3utr]
        coding_region_seqs.append(coding_region)
            
    # Write the sequences to separate files
    SeqIO.write(five_utr_seqs, args.output_5UTR, "fasta")
    SeqIO.write(coding_region_seqs, args.output_coding_region, "fasta")
    SeqIO.write(three_utr_seqs, args.output_3UTR, "fasta")
            