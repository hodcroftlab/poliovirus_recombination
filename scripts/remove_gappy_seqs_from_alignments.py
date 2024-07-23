### This script removes sequences with more than the given proportion of gaps/Ns from the alignment.

from Bio import SeqIO
import argparse

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Remove sequences with a proportion of gaps higher than max_gap_fraction from the alignment.")
    argparser.add_argument("--alignments", help="Input alignment file(s) in fasta format.", nargs="+")
    argparser.add_argument("--max_gap_fraction", help="Maximum fraction of gaps allowed in a sequence.", type=float)
    argparser.add_argument("--output_dir", help="Output directory for the filtered alignments.")
    
    args = argparser.parse_args()
    
    # Read the alignment
    alignments = {}
    for alignment_file in args.alignments:
        name = alignment_file.split("/")[-1].split("_")[0]
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        # Add the alignment to the dictionary
        alignments[name] = alignment

    # Filter alignments based on the maximum gap fraction
    filtered_alignments = {}
    for name, alignment in alignments.items():
        filtered_alignment = [rec for rec in alignment if sum(char == "-" for char in rec.seq) / len(rec.seq) <= args.max_gap_fraction]
        filtered_alignments[name] = filtered_alignment
    
    # Write the filtered alignments to the output directory
    for name, alignment in filtered_alignments.items():
        output_file = f"{args.output_dir}/{name}_filtered.fasta"
        SeqIO.write(alignment, output_file, "fasta")