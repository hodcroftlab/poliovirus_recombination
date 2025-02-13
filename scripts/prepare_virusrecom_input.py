### This script prepares the input files for running VirusRecom.

# Required input files for VirusRecom:
# 1. Alignment file (.fasta); each sequence name requires a lineage/serotype mark
# 2. Text file of lineage/serotype marks (.txt); simple listing, one per line

import pandas as pd
from Bio import SeqIO
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare input files for VirusRecom.")
    parser.add_argument("--alignment", help="Path to fasta file containing sequences.")
    parser.add_argument("--groups", nargs="+", help="Groups to consider recombiants in.")
    parser.add_argument("--metadata", help="Path to metadata file (csv/tsv).")
    parser.add_argument("--metadata_id_col", help="Column name in metadata file containing sequence identifiers matching the fasta files (e.g., accession IDs).")
    parser.add_argument("--metadata_group_col", help="Column name in metadata file containing assignment of sequences to reference groups (e.g., serotype, clade membership).")
    parser.add_argument("--output_alignment", help="Path to output alignment file (.fasta).")
    parser.add_argument("--output_reference_txt", help="Path to output .txt file listing the reference groups.")
    
    args = parser.parse_args()

    # Read alignment
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    # alignment = list(SeqIO.parse("../data/sequences/alignments/filtered/fullgenome_filtered.fasta", "fasta"))
    
    # Read metadata for serotype information
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
        
    # metadata = pd.read_csv("../data/metadata/evc_full_genomes_metadata_rivm.csv")
    
    # For Sabin sequences, change the serotype to "OPV1/2/3"
    sabin_dict = {"AY184219.1": "OPV1", "AY184220.1": "OPV2", "AY184221.1": "OPV3"} 
    for seq_id, serotype in sabin_dict.items():
        metadata.loc[metadata[args.metadata_id_col] == seq_id, args.metadata_group_col] = serotype
    
    
    # Create new alignment with sequences named as serotype mark + accession ID
    new_alignment = []
    
    for seq in alignment:
        # Check if the sequence is in the metadata file
        if seq.id not in metadata[args.metadata_id_col].values:
            print(f"Sequence {seq.id} not found in metadata file column '{args.metadata_id_col}'.")
            continue
        
        serotype_mark = metadata.loc[metadata[args.metadata_id_col] == seq.id, args.metadata_group_col].values[0]  # Get serotype mark for sequence
        seq.id = serotype_mark + "_" + seq.id  # Update sequence ID
        new_alignment.append(seq)
        new_alignment[-1].description = ""
    
    # Write new alignment
    SeqIO.write(new_alignment, args.output_alignment, "fasta")


    # Write .txt file listing the unique group assignments (e.g., clade_membership), one per line
    with open(args.output_reference_txt, "w") as f:
        for group in metadata[args.metadata_group_col].unique():
            f.write(group + "\n")
    
    
    