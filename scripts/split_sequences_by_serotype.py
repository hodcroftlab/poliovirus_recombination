### This script takes alignment as well as metadata files as inputs and splits the sequences in the alignment file into clusters based on the metadata.

from Bio import SeqIO
import pandas as pd
import argparse

if __name__ == "__main__":
    
    # Parse command-line arguments
    argparser = argparse.ArgumentParser(description="Split sequences in an alignment file into clusters based on metadata.")
    argparser.add_argument("--alignments", nargs="+", help="Input alignment file(s) in fasta format.")
    argparser.add_argument("--metadata", help="Metadata file in CSV/TSV format.")
    argparser.add_argument("--metadata_id_col", help="Column in the metadata file containing sequence identifiers (accession IDs).")
    argparser.add_argument("--metadata_type_col", help="Column in the metadata file to use for grouping.")
    argparser.add_argument("--outdir", help="Output directory for the split alignment files.")
    args = argparser.parse_args()
    
    # Read the metadata file
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
    else:
        raise ValueError("Metadata file must be a CSV or TSV file.")
    
    # Read the alignment file(s)
    for alignment_file in args.alignments:
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        region = alignment_file.split("/")[-1].split("_")[0]
        
        print(f"\nProcessing {region} alignment ...")
        
        # Group sequences by serotype
        for serotype in metadata[args.metadata_type_col].unique():
            serotype_name = serotype.replace("/", "_")
            print(f"Cluster: {serotype_name}")
            serotype_alignment = [seq for seq in alignment if seq.id in metadata[metadata[args.metadata_type_col] == serotype][args.metadata_id_col].values]
            
            # Write the serotype alignment to a new file
            output_file = f"{args.outdir}/{serotype_name}/{serotype_name}_{region}_aligned.fasta"
            SeqIO.write(serotype_alignment, output_file, "fasta")
        