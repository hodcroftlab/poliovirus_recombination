### This script calculates p-distances between chunks of each sequence and the consensus sequences for the serotypes that are given.

import pandas as pd
import argparse
import numpy as np
from Bio import SeqIO

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Calculate p-distances between chunks of each sequence and the consensus sequences for the serotypes that are given.")
    argparser.add_argument("--alignment", help="Input alignment file in fasta format.")
    argparser.add_argument("--consensus", help="Input consensus sequence file in fasta format (can contain multiple sequences).")
    argparser.add_argument("--metadata", help="Metadata file (.csv/.tsv) containing sequence identifiers and group information.")
    argparser.add_argument("--metadata_id_col", help="Column name in the metadata file containing the sequence identifiers (default: strain).", default="strain")
    argparser.add_argument("--metadata_group_col", help="Column name in the metadata file containing the group identifiers (default: serotype).", default="serotype")
    argparser.add_argument("--chunk_size", type=int, help="Size of the chunks to compare with the consensus sequences.")
    argparser.add_argument("--output", help="Output file (.csv/.tsv) for a dataframe containing the calculated pairwise distances.")
    args = argparser.parse_args()
    
    # Load data files
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    consensus = list(SeqIO.parse(args.consensus, "fasta"))
    
    # alignment = list(SeqIO.parse("../data/sequences/alignments/full_genome/assembled/fullgenome_pre_masking_assembled.fasta", "fasta"))
    # consensus = list(SeqIO.parse("../data/sequences/consensus_sequences/EVC109_consensus.fasta", "fasta"))
    # chunk_size = 200
    
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
        
    metadata = metadata[[args.metadata_id_col, args.metadata_group_col]]
    
    # Parameters
    chunk_size = args.chunk_size
    metadata_id_col = args.metadata_id_col
    metadata_group_col = args.metadata_group_col
    
    # Split alignment and consensus sequences into chunks
    seq_len_alignment = len(alignment[0].seq)
    seq_len_consensus = len(consensus[0].seq)
    if seq_len_alignment != seq_len_consensus:
        raise ValueError("The alignment and consensus sequences have different lengths.")
    
    chunks = {}
    for i in range(0, seq_len_alignment, chunk_size):
        chunk_alignment = []
        for record in alignment:
            chunk_alignment.append(record[i:i+chunk_size])
        chunk_consensus = []
        for record in consensus:
            chunk_consensus.append(record[i:i+chunk_size])
        
        # Adjust the chunk size for the last chunk
        if i + chunk_size > seq_len_alignment:
            chunk_size = seq_len_alignment - i
        name = f"{i+1}-{i+chunk_size}"
        chunks[name] = (chunk_alignment, chunk_consensus) # Store the chunks in a dictionary
     
    # For each sequence in the alignment, calculate pairwise distance to each sequence in the consensus alignment
    # Initialize empty list to collect results
    results = []
    for chunk, (alignment_chunk, consensus_chunk) in chunks.items():
        print(f"Calculating p-distances for chunk {chunk} ...")
        for record in alignment_chunk:
            seq1 = record.seq 
            for consensus_record in consensus_chunk:
                
                seq2 = consensus_record.seq
                
                # Convert to numpy arrays for faster element-wise calculations
                seq1 = np.array(list(seq1))
                seq2 = np.array(list(seq2))
                
                # Mask all positions where either of the two sequences has a gap or an N
                valid_positions = (seq1 != "-") & (seq1 != "N") & (seq2 != "-") & (seq2 != "N")
                seq1_valid = seq1[valid_positions]
                seq2_valid = seq2[valid_positions]
                
                # Update the sequence length
                seq_len_valid = len(seq1_valid)
                # if seq_len_valid < (chunk_size/10):
                #     raise Warning(f"For {record.id} and {consensus_record.id}, the number of valid positions is less than 10% of the chunk size.")
                
                # Now calculate the number of differing positions
                nd = np.sum(seq1_valid != seq2_valid)
                
                # Calculate the p-distance
                
                if seq_len_valid == 0:
                    pdist = np.nan
                else:
                    pdist = round(nd / seq_len_valid, 4)
                
                # Append the result as a tuple to the results list
                results.append((chunk, record.id, consensus_record.id, pdist))
                
    results_df = pd.DataFrame(results, columns=["chunk", "query_seq", "consensus_seq", "pdist"]) # Convert the list of results to a dataframe
    results_df = results_df.merge(metadata, left_on="query_seq", right_on=metadata_id_col) # Merge with metadata to get serotype information for the query_seqs
    results_df.rename(columns={metadata_group_col: "assigned_serotype"}, inplace=True) # Rename metadata_group_col to assigned_serotype
    results_df = results_df[["chunk", "query_seq", "assigned_serotype", "consensus_seq", "pdist"]] # Reorder columns
    results_df["start"] = results_df["chunk"].apply(lambda x: int(x.split("-")[0])) # Split the "chunk" column into numeric "start" value for correct sorting
    results_df.sort_values(by=["assigned_serotype", "query_seq", "start", "pdist"], inplace=True) # Order rows by assigned_serotype, query_seq, start, and pdist
    results_df.drop(columns=["start"], inplace=True) # Drop the "start" column
    
    # Additionally, output filtered dataframe that only contains the lowest three p-distances for each chunk-query_seq pair
    results_df_filtered = results_df.groupby(["chunk", "query_seq"]).apply(lambda x: x.nsmallest(3, "pdist")).reset_index(drop=True)
    results_df_filtered["start"] = results_df_filtered["chunk"].apply(lambda x: int(x.split("-")[0])) # Split the "chunk" column into numeric "start" value for correct sorting
    results_df_filtered.sort_values(by=["assigned_serotype", "query_seq", "start", "pdist"], inplace=True) # Order rows by assigned_serotype, query_seq, start, and pdist
    results_df_filtered.drop(columns=["start"], inplace=True) # Drop the "start" column
    
    # Save the results to a file
    if args.output.endswith(".csv"):
        results_df.to_csv(args.output, index=False)
        results_df_filtered.to_csv(args.output.replace(".csv", "_filtered.csv"), index=False)
    elif args.output.endswith(".tsv"):
        results_df.to_csv(args.output, index=False, sep="\t")
        results_df_filtered.to_csv(args.output.replace(".tsv", "_filtered.tsv"), index=False)
    
    