### This script filters the dataframe of neighbors for reduce the number of sequences for the next steps.

# Filtering strategy:
# 1. Filter out any neighbors that have the same serotype as the query sequence
# 2. Remove sequence-region pairs where the closest neighbor has a distance greater than a given distance_threshold
# This is primarily to reduce the number of sequences to a subset with higher confidence of finding recombination partners

import pandas as pd
import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Filter the dataframe of neighbors.")
    argparser.add_argument("--neighbors", help="Input file (CSV/TSV) containing the dataframe of neighbors.")
    argparser.add_argument("--distance_threshold", type=float, help="Distance threshold to filter out sequence-region pairs.", default=0.1)
    argparser.add_argument("--out", help="Output file (CSV/TSV) containing the filtered dataframe of neighbors.")
    args = argparser.parse_args()


    # Read the dataframe of neighbors
    if args.neighbors.endswith(".csv"):
        neighbors = pd.read_csv(args.neighbors)
    elif args.neighbors.endswith(".tsv"):
        neighbors = pd.read_csv(args.neighbors, sep="\t")
    else:
        raise ValueError("Neighbors file must be in CSV or TSV format.")
    
    #neighbors = pd.read_csv("../data/recombination_tools/custom_simplots/neighbors.csv")

    #-------------------------------------------------------------------------------------------------------------------------------------------------
    # Step 1: Filter out any neighbors that have the same serotype as the query sequence
    #-------------------------------------------------------------------------------------------------------------------------------------------------

    print("Step 1: Filtering out neighbors that have the same serotype as the query sequence...")
    neighbors_filtered = neighbors[neighbors["serotype"] != neighbors["neighbor_serotype"]]

    #-------------------------------------------------------------------------------------------------------------------------------------------------
    # Step 2: Filter out neighbors that have a distance greater than a given distance_threshold
    #-------------------------------------------------------------------------------------------------------------------------------------------------

    distance_threshold = args.distance_threshold
    #distance_threshold = 0.1

    print(f"Step 2: Filtering out neighbors that have a distance greater than {distance_threshold}...")
    neighbors_filtered = neighbors_filtered[neighbors_filtered["distance"] <= distance_threshold]
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------
    # Step 3: For each sequence-region pair, keep only the closest neighbor
    #-------------------------------------------------------------------------------------------------------------------------------------------------

    print("Step 3: Keeping only the closest neighbor for each sequence-region pair...")
    
    # Create new column concatenating seq and region
    neighbors_filtered["seq_region"] = neighbors_filtered["seq"] + "_" + neighbors_filtered["region"]

    # Keep the closest neighbor for each sequence-region pair
    neighbors_filtered = neighbors_filtered[neighbors_filtered["distance"] == neighbors_filtered.groupby("seq_region")["distance"].transform("min")]
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------
    # Step 4: If a sequences has two closest neighbors with the same distance, keep only one of them
    #-------------------------------------------------------------------------------------------------------------------------------------------------

    print("Step 4: Removing duplicates (if any)...")
    neighbors_filtered = neighbors_filtered.drop_duplicates(subset="seq_region")
    
    # Remove seq_region column
    neighbors_filtered.drop(columns=["seq_region"], inplace=True)
    
    print(f"Number of rows before filtering: {neighbors.shape[0]}")
    print(f"Number of rows after filtering: {neighbors_filtered.shape[0]}")

    # Save the filtered neighbors dataframe
    if args.out.endswith(".csv"):
        neighbors_filtered.to_csv(args.out, index=False)
    elif args.out.endswith(".tsv"):
        neighbors_filtered.to_csv(args.out, sep="\t", index=False)
    else:  
        raise ValueError("Output file must be in CSV or TSV format.")