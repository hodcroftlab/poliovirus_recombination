### This script compares p-distances between VP1 and the specified regions and applies a cutoff to the absolute difference in p-distances to identify potential recombinant sequences.

import pandas as pd
import argparse
import numpy as np

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Identify potential recombinant sequences based on p-distance comparisons.")
    argparser.add_argument("--pdistances", nargs="+", help="File(s) containing pairwise p-distances.")
    argparser.add_argument("--regions", nargs="+", help="Regions to compare with VP1.")
    argparser.add_argument("--cutoff", type=float, help="Cutoff value for p-distance difference between VP1 and the selected regions to identify region as potentially recombinant.")
    argparser.add_argument("--output", help="Output file (.csv) for the candidate recombinant sequences.")
    argparser.add_argument("--plotdir", help="Output directory for histograms of sequence occurences beyond the cutoff.")
    args = argparser.parse_args()
    
    # Load data files
    pdistances = args.pdistances
    regions = args.regions
    cutoff_candidates = args.cutoff    
    
    pdistance_dict = {} # Prepare the dictionary with dataframes
    for file in pdistances:
        serotype = file.split("/")[-1].replace("_pdistances.csv", "")
        pdistance_dict[serotype] = pd.read_csv(file)
    
    # Initialize results dataframes
    results = pd.DataFrame(columns=["seq1", "seq2", "diff", "serotype", "region"])
        
    # For each serotype, for each region, compare p-distances with VP1
    for serotype, df in pdistance_dict.items():
        for region in regions:
            # Subset dataframe for VP1 and the specified region
            df_subset = df[["seq1", "seq2", "VP1", region]].copy()
            
            # Calculate the absolute difference in p-distances
            # Adjust calculation for 5'UTR --> assumption: more conserved than other regions
            if region == "5UTR":
                df_subset["diff"] = abs(df_subset["VP1"]*0.8 - df_subset[region])
            else:
                df_subset["diff"] = abs(df_subset["VP1"] - df_subset[region])
            
            # Filter for sequences with a differences above/below the cutoffs
            df_candidates = df_subset[df_subset["diff"] > cutoff_candidates].copy()
            
            # Add serotype and region columns for subsequent merging
            df_candidates["serotype"] = serotype
            df_candidates["region"] = region
            
            # Remove VP1 and region columns to avoid redundancy
            df_candidates.drop(columns=["VP1", region], inplace=True)
            
            # Add to results dataframe
            results = pd.concat([results, df_candidates])
                        
    # Concatenate seq1 and seq2 into a single column for counting
    results_long = pd.melt(results, id_vars=["serotype", "region"], value_vars=["seq1", "seq2"], var_name="seq_type", value_name="seq")

    # Group by serotype, region, and sequence, then count occurrences
    counts = results_long.groupby(["serotype", "region", "seq"]).size().reset_index(name="count")
    
    # If total sequences in a serotype and region is <5, remove this serotype and region from the dataframe
    counts = counts[counts.groupby(["serotype", "region"])["count"].transform("sum") >= 5]
    
    # Create a new column: divide the count for each sequence (within a serotype and region) by the total number of sequences within that serotype and region
    def calculate_proportion(row):
        serotype = row["serotype"]
        region = row["region"]
        # Count rows where the serotype and region match
        total_seqs = counts[(counts["serotype"] == serotype) & (counts["region"] == region)].shape[0]
        if total_seqs == 0:
            return np.nan
        else:
            proportion = row["count"] / total_seqs
            return round(proportion, 3)
    
    counts["proportion_matched"] = counts.apply(calculate_proportion, axis=1)
    
    # Filter for sequences that match more than 50% of the sequences that are beyond the cutoff
    candidate_recombinants = counts[counts["proportion_matched"] >= 0.5]
    
    # For seqs with occurences in multiple regions, combine the regions into a single string
    
    # Convert "count" and "proportion_matched" columns to string before aggregation
    candidate_recombinants = candidate_recombinants.copy()
    candidate_recombinants[["count", "proportion_matched"]] = candidate_recombinants[["count", "proportion_matched"]].astype(str)

    # Group by and aggregate
    candidate_recombinants = candidate_recombinants.groupby(["seq", "serotype"]).agg({
        "region": ", ".join,
        "count": ", ".join,
        "proportion_matched": ", ".join
    }).reset_index()
    
    
    # Sort by serotype and seq
    candidate_recombinants = candidate_recombinants.sort_values(by=["serotype", "seq"])
    
    # Save the results and counts dataframes to file
    candidate_recombinants.to_csv(args.output, index=False)
    results.to_csv(args.output.replace(".csv", "_diff.csv"), index=False)
    counts.to_csv(args.output.replace(".csv", "_counts.csv"), index=False)    
    
    # Plot distribution of sequence occurences above the threshold per serotype and region
    # For each serotype, make a grid of histograms (e.g., 2x2 for 4 regions)
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    
    regions.sort()
    n_regions = len(regions)
    
    # Check if plotdir exists, if not create it
    import os
    if not os.path.exists(args.plotdir):
        os.makedirs(args.plotdir)
    
    for serotype in counts["serotype"].unique():
        fig, axs = plt.subplots(n_regions // 2, 2, figsize=(10, 10))
        for i, region in enumerate(regions):
            ax = axs[i // 2, i % 2]
            data = counts[(counts["serotype"] == serotype) & (counts["region"] == region)].copy()
            
            # Ensure there is data for the histogram
            if not data.empty:
                # Calculate min and max for binning
                min_data = np.floor(data["count"].min()) - 0.5
                max_data = np.ceil(data["count"].max()) + 0.5
                
                # Calculate the range of the data
                data_range = max_data - min_data
                
                # Set base bin width to 1 if data range is small, otherwise increase it
                if data_range < 50:
                    bin_width = 1  # Default bin width for small ranges
                else:
                    # Increase the bin width based on the data range for large datasets
                    bin_width = data_range // 30  # Aim for ~30 bins

                # Check if min_data is less than max_data before creating bins
                if min_data < max_data:
                    # Create bins centered around integers
                    bins = np.arange(min_data, max_data + bin_width, bin_width)

                    # Plot histogram with the new bins
                    sns.histplot(data["count"], bins=bins, ax=ax)
                else:
                    sns.histplot(data["count"], bins=1, ax=ax)
                    
                # Add horizontal line at 50% of the total unique sequences above the threshold
                total_seqs = data.shape[0]
                threshold = total_seqs / 2
                ax.axvline(x=threshold, color="red", linestyle="--")
                
                # Set x- and y-axis to have ticks only at whole numbers
                ax.xaxis.set_major_locator(MaxNLocator(integer=True)) 
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                
            else:
                # Handle empty data case: Set limits and ticks to show only 0
                ax.set_xlim([-0.5, 0.5])
                ax.set_ylim([-0.5, 0.5])  # Minimal height for the y-axis
                ax.set_xticks([0])
                ax.set_yticks([0])
                ax.set_xlabel("Sequence Matches")
                ax.set_ylabel("Count")    
                
            serotype_name = serotype.replace("_", "/")
            ax.set_title(f"{serotype_name} - {region}")
            ax.set_xlabel("Sequence Matches")
            ax.set_ylabel("Count")

            
        plt.tight_layout()
        plt.savefig(f"{args.plotdir}/{serotype}_sequence_counts_histograms.png")
        plt.close()
        
    
    
    