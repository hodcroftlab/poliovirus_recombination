### This script takes phylogenetic trees in newick format as input and finds the nearest neighbors of candidate recombinants.
### Search the trees for the regions that are alleged recombinants.

from Bio import Phylo
import argparse
import pandas as pd

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def find_neighbors_for_candidates(candidates, tree_dict, metadata, n_neighbors):
    """
    For each candidate recombinant, find the nearest neighbors by moving upwards from the terminal node
    to common ancestors, stopping when a clade containing at least n_neighbors terminal nodes is found.
    """
    results = []  # Use a list to collect results

    for _, row in candidates.iterrows():
        seq = row["seq"]
        region = row["region"]
        serotype = row["serotype"]

        if region not in tree_dict:
            print(f"Warning: Tree for region {region} not found. Skipping.")
            continue

        tree = tree_dict[region]

        # Find the terminal node (candidate sequence)
        candidate_node = None
        for terminal in tree.get_terminals():
            if terminal.name == seq:
                candidate_node = terminal
                break

        if candidate_node is None:
            # Skip the sequence for this region
            print(f"Sequence {seq} not found in tree for region {region}. Skipping.")
            continue

        # Traverse upwards to find a clade with enough terminals
        print(f"Finding neighbors for {seq} in region {region}...")
        current_clade = candidate_node
        while current_clade:
            terminals_in_clade = current_clade.get_terminals()

            if len(terminals_in_clade) >= n_neighbors + 1:  # +1 to exclude the candidate itself
                break

            # Move upwards to the parent
            current_clade = get_parent(tree, current_clade)
            

        # Now we have enough terminals; compute distances within this clade
        print("Computing distances...")
        neighbors = []
        for terminal in terminals_in_clade:
            if terminal != candidate_node:
                distance = tree.distance(candidate_node, terminal)
                neighbors.append((terminal.name, distance))

        # Sort neighbors by distance and select the closest ones
        neighbors = sorted(neighbors, key=lambda x: x[1])[:n_neighbors]

        # Append results to the list
        for neighbor, distance in neighbors:
            results.append([seq, serotype, region, neighbor, distance])

    # Convert the results list to a DataFrame at the end
    output_df = pd.DataFrame(results, columns=["seq", "serotype", "region", "neighbor_accn", "distance"])
    
    # Merge with metadata to get serotype information
    output_df = output_df.merge(metadata, left_on="neighbor_accn", right_on=args.metadata_id_col, how="left").drop(columns=[args.metadata_id_col]).rename(columns={args.metadata_type_col: "neighbor_serotype"})
    
    # Sort by seq, region, and distance
    output_df = output_df.sort_values(by=["serotype", "seq", "region", "distance"])

    return output_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find nearest neighbors of candidate recombinants in phylogenetic trees.")
    parser.add_argument("--trees", nargs="+", help="Input Newick files containing phylogenetic trees.", required=True)
    parser.add_argument("--candidates", help="CSV/TSV file of candidate recombinants and regions.")
    parser.add_argument("--metadata", help="Metadata file (CSV/TSV) containing serotype information for all sequences.")
    parser.add_argument("--metadata_id_col", help="Column name in the metadata file containing sequence identifiers (accession IDs).")
    parser.add_argument("--metadata_type_col", help="Column name in the metadata file containing serotypes.")
    parser.add_argument("--n_neighbors", type=int, help="Number of nearest neighbors to find.", default=3)
    parser.add_argument("--out", help="Output file (CSV/TSV) containing nearest neighbors.", required=True)
    args = parser.parse_args()
    
    # Read metadata
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
    else:
        raise ValueError("Metadata file must be in CSV or TSV format.")
    
    # Select columns we need
    metadata = metadata[[args.metadata_id_col, args.metadata_type_col]]
    
    #metadata = pd.read_csv("../data/metadata/evc_full_genomes_metadata_rivm.csv")
    #metadata = metadata[["Accession", "serotype"]]

    # If candidate argument is provided, read candidates and find neighbors
    if args.candidates:
        print("Using provided candidates.\n")
        # Read candidate recombinants
        if args.candidates.endswith(".csv"):
            candidates = pd.read_csv(args.candidates)
        elif args.candidates.endswith(".tsv"):
            candidates = pd.read_csv(args.candidates, sep="\t")
        else:
            raise ValueError("Candidates file must be in CSV or TSV format.")
        
        # Select only the columns we need
        candidates = candidates[["seq", "serotype", "region"]]

        # Split regions into separate rows
        candidates = candidates.assign(region=candidates["region"].str.split(", ")).explode("region")
        
        # Add VP1 as an extra row for every unique sequence
        vp1_candidates = candidates.groupby("seq").first().reset_index()
        vp1_candidates["region"] = "VP1"
        candidates = pd.concat([candidates, vp1_candidates], ignore_index=True)
        
    else:
        print("No candidates provided. Using all sequences in the metadata as candidates.\n")
        # Use all sequences as candidates; get accessions from metadata
        # Use pre-defined regions: VP1, 2C, 3D, 5UTR
        candidates = metadata.copy()
        candidates = candidates.rename(columns={args.metadata_id_col: "seq", args.metadata_type_col: "serotype"}) # Rename "Accession" to "seq"
        # Elongate the dataframe by adding all regions to each sequence
        regions = ["VP1", "2C", "3D", "5UTR"]
        
        # Repeat each row 4 times
        candidates = candidates.loc[candidates.index.repeat(len(regions))].reset_index(drop=True)

        # Assign genomic regions to the new column
        candidates["region"] = regions * len(metadata)


    # Read phylogenetic trees into a dictionary based on region
    tree_dict = {tree.split("/")[-1].split("_")[0]: Phylo.read(tree, "newick") for tree in args.trees}
    

    # Find neighbors for all candidates
    output = find_neighbors_for_candidates(candidates, tree_dict, metadata, args.n_neighbors)

    # Write output to file
    if args.out.endswith(".csv"):
        output.to_csv(args.out, index=False)
    elif args.out.endswith(".tsv"):
        output.to_csv(args.out, sep="\t", index=False)
    else:
        raise ValueError("Output file must be in CSV or TSV format.")