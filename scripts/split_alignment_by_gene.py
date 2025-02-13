### This script will take an aligned fasta file and split it into multiple files containing the sequences of the individual genes; based on gene annotations in the gb file

from Bio import SeqIO
import pandas as pd
import argparse


if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description="Split a multiple sequence alignment (.fasta) into several alignments for each gene based on gene annotations in a .csv file.")
    argparser.add_argument("--alignment", help="Input aligned fasta file to be split into .")
    argparser.add_argument("--gene_annotation", help="Gene annotation file in CSV or TSV format (expected columns: feature, start, end).")
    argparser.add_argument("--output_dir", help="Output directory for the split alignments.")
    
    args = argparser.parse_args()

    # Read in aligned fasta
    alignment = list(SeqIO.parse(args.alignment, "fasta"))

    # Read in gene annotations
    if args.gene_annotation.endswith(".csv"):
        gene_annotations = pd.read_csv(args.gene_annotation)
    elif args.gene_annotation.endswith(".tsv"):
        gene_annotations = pd.read_csv(args.gene_annotation, sep="\t")
    else:
        print("Gene annotation file must be in CSV or TSV format.")

    # Split aligned fasta by gene
    # Iterate over genes/features
    for feature in gene_annotations["feature"]:
        gene_alignment = []
        start = gene_annotations[gene_annotations["feature"] == feature]["start"].values[0]-1
        end = gene_annotations[gene_annotations["feature"] == feature]["end"].values[0]
        # Iterate over records in aligned fasta
        for record in alignment:
            gene_record = record[start:end]
            gene_alignment.append(gene_record)
        # Write gene records to new fasta file
        SeqIO.write(gene_alignment, f"{args.output_dir}/{feature}_aligned.fasta", "fasta")
            