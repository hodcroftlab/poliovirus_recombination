### This script extracts sub-alignments from the input alignment for each group of sequences specified in the metadata.

import argparse
import pandas as pd
from Bio import SeqIO

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Extract sub-alignments from the input alignment for each group of sequences specified in the metadata.")
    argparser.add_argument("--alignment", help="Input file containing aligned sequences.")
    argparser.add_argument("--metadata", help="CSV/TSV file containing metadata for the aligned sequences.")
    argparser.add_argument("--metadata_id_col", help="Column name in the metadata file containing the sequence identifiers.")
    argparser.add_argument("--metadata_group_col", help="Column name in the metadata file containing the group identifiers.")
    argparser.add_argument("--outdir", help="Output directory to write the sub-alignments.")
    args = argparser.parse_args()
    
    # Read the alignment
    alignment = SeqIO.parse(args.alignment, "fasta")
    alignment_dict = {record.id: record for record in alignment}
    
    # Read the metadata
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep = "\t")
    else:
        raise ValueError("Metadata file must be a CSV or TSV file.")
    
    # Group the sequences by the group column
    group_dict = metadata.groupby(args.metadata_group_col)[args.metadata_id_col].apply(list).to_dict()
    
    # Extract sub-alignments for each group
    for group, group_sequences in group_dict.items():
        group = group.replace("/", "_")
        sub_alignment = [alignment_dict[seq_id] for seq_id in group_sequences]
        output_file = f"{args.outdir}/{group}.fasta"
        SeqIO.write(sub_alignment, output_file, "fasta")
        
    # In addition, write all and non-polio sequences to separate files
    non_polio_groups = [cluster for cluster in metadata[args.metadata_group_col].unique().tolist() if cluster not in ["PV1", "PV2", "PV3"]]
    non_polio_alignment = [alignment_dict[seq_id] for seq_id in metadata[metadata[args.metadata_group_col].isin(non_polio_groups)][args.metadata_id_col]]
    SeqIO.write(non_polio_alignment, f"{args.outdir}/non_polio.fasta", "fasta")
    
    SeqIO.write(list(alignment_dict.values()), f"{args.outdir}/all.fasta", "fasta")