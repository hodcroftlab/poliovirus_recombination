### This script creates sub-alignments for each candidate recombinant sequence, its neighbors, and the consensus sequences.
import pandas as pd
import argparse
from Bio import SeqIO
import os

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Create alignments for custom SimPlots.")
    argparser.add_argument("--alignment", help="Alignment file (fasta).")
    argparser.add_argument("--neighbors", help="CSV/TSV containing the neighbors.")
    argparser.add_argument("--neighbors_candidate_col", help="Column in the neighbors CSV/TSV containing candidate recombinants.")
    argparser.add_argument("--neighbors_neighbor_col", help="Column in the neighbors CSV/TSV containing neighbor accessions.")
    argparser.add_argument("--consensus", help="Consensus sequences (fasta).")
    argparser.add_argument("--sabin", nargs="+", help="Accession numbers of Sabin sequences (1, 2, and 3) in the input alignment.")
    argparser.add_argument("--outdir", help="Output directory.")
    args = argparser.parse_args()

    # Read alignment
    alignment = list(SeqIO.parse(args.alignment, "fasta"))

    # Read filtered neighbors
    if args.neighbors.endswith(".csv"):
        neighbors = pd.read_csv(args.neighbors)
    elif args.neighbors.endswith(".tsv"):
        neighbors = pd.read_csv(args.neighbors, sep="\t")
    else:
        raise ValueError("Please provide a CSV or TSV file for the neighbors.")

    # For each seq in neighbors, get all the neighbor accesions and save in dictionary
    neighbor_dict = {}
    for seq in neighbors[args.neighbors_candidate_col].unique():
        neighbors_filtered = neighbors[neighbors[args.neighbors_candidate_col] == seq]
        neighbor_accn = neighbors_filtered[args.neighbors_neighbor_col].unique()
        neighbor_dict[seq] = neighbor_accn
        
    # Read consensus sequences
    consensus = list(SeqIO.parse(args.consensus, "fasta"))
    
    # Get the Sabin sequences
    sabin_accessions = args.sabin
    sabin_sequences = [record for record in alignment if record.id in sabin_accessions]
        
    # Build an alignment for each seq, its neighbors, the consensus sequences, and the Sabin sequences
    alignment_dict = {}
    for seq, neighbor_accn in neighbor_dict.items():
        alignment_dict[seq] = [record for record in alignment if record.id in [seq] + list(neighbor_accn)] # Add candidate and neighbor sequences
        alignment_dict[seq] += [record for record in consensus] # Add consensus sequences
        alignment_dict[seq] += sabin_sequences # Add Sabin sequences
        
        # Place the candidate sequence at the start
        alignment_dict[seq] = [record for record in alignment_dict[seq] if record.id == seq] + [record for record in alignment_dict[seq] if record.id != seq]

    # Check if outdir exists, if not create it
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Save alignments
    for seq, records in alignment_dict.items():
        SeqIO.write(records, f"{args.outdir}/{seq}.fasta", "fasta")
    


