### This script counts the number of sequences in a fasta file that have the same length and reports the counts for each length.

from Bio import SeqIO
import argparse

argparser = argparse.ArgumentParser(description="Count the number of sequences in a fasta file that have the same length and report the counts for each length.")
argparser.add_argument("--fasta", help="Input fasta file")
args = argparser.parse_args()

lengths = {}
for record in SeqIO.parse(args.fasta, "fasta"):
    length = len(record.seq)
    if length in lengths:
        lengths[length] += 1
    else:
        lengths[length] = 1
        
for length in sorted(lengths.keys()):
    print("%d\t%d" % (length, lengths[length]))