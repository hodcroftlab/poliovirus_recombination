### Plot pdistance of each genomic region against other genomic regions in a dataframe (contains one colummn for each genomic region)

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
import numpy as np
import seaborn as sns
import pandas as pd
import re
import os
import argparse
import mpl_scatter_density # adds projection='scatter_density'

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Plot pairwise p-distances between VP1 and 3Dpol.")
    argparser.add_argument("--pdistances", nargs="+", help="Input file(s) (.csv) containing pairwise p-distances.")
    argparser.add_argument("--reference_region", help="Genomic region to use as reference (plotted on the x-axis).")
    argparser.add_argument("--compare_regions", nargs="+", help="Genomic regions to compare against the reference region.")
    argparser.add_argument("--serotypes", nargs="+", help="Serotypes to compare.")
    argparser.add_argument("--grid", action="store_true", help="Make grid of subplots (serotypes x region comparisons).")
    argparser.add_argument("--outdir", help="Output directory for the plots.")
    argparser.add_argument("--outformat", help = "Output format for the plots (options: png/pdf/svg).")
    
    args = argparser.parse_args()
    
    # Set font to Arial
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 30
    
    # Delete the output directory if it already exists
    outdir = args.outdir
    if os.path.exists(outdir):
        os.system(f"rm -r {outdir}")
    os.makedirs(outdir)
    
    
    serotypes = args.serotypes
    reference_region = args.reference_region
    compare_regions = args.compare_regions
    
    print(f"Serotypes: {serotypes}")
    print(f"Reference region: {reference_region}")
    print(f"Regions to compare: {compare_regions}")
    
    
    # Load data files
    pdistances = args.pdistances
    pdistance_dict = {} # Prepare the dictionary with dataframes
    for file in pdistances:
        serotype = file.split("/")[-1].replace("_pdistances.csv", "")
        if serotype in serotypes:
            pdistance_dict[serotype] = pd.read_csv(file)
            
    # Exclude weird PV1 sequences until further investigation
    pv1_exclude = ["OR208593.1", "OR208600.1", "OR208605.1"]
    for serotype, df in pdistance_dict.items():
        df = df[~df["seq1"].isin(pv1_exclude)]
        df = df[~df["seq2"].isin(pv1_exclude)]
        pdistance_dict[serotype] = df
    
    # Find maximum difference between p-distances for each serotype to set the color scale
    max_diff = []
    for serotype, df in pdistance_dict.items():
        for region in compare_regions:
            df_subset = df[["seq1", "seq2", reference_region, region]].copy()
            df_subset = df_subset.dropna()
            df_subset["diff"] = abs(df_subset[region] - df_subset[reference_region])
            max_diff.append(df_subset["diff"].max())
    max_diff = max(max_diff)
    print("Generating plots ...")
    
    # Define the color normalization for consistent color mapping
    norm = mcolors.Normalize(vmin=0, vmax=max_diff)  # Set the range of hue values
    cmap = plt.get_cmap("viridis")  # Color palette for the scatterplots
    
    
    if args.grid:
        # Make grid of subplots (serotypes x region comparisons)
        n_serotypes = len(pdistance_dict)
        n_regions = len(compare_regions)
        fig, axes = plt.subplots(n_regions, n_serotypes, figsize=(8*n_serotypes, 8*n_regions))
        
        letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")  # Unique letters for subplots
        letter_count = 0  # Counter for letters

        # Loop over each subplot in the grid
        for i, serotype in enumerate(pdistance_dict.keys()):
            for j, region in enumerate(compare_regions):
                ax = axes[j, i]  # Get the axis for the current serotype and region
                    
                # Scatter plot
                df_subset = pdistance_dict[serotype][["seq1", "seq2", reference_region, region]].copy()
                df_subset = df_subset.dropna()
                df_subset["diff"] = abs(df_subset[region] - df_subset[reference_region])
                
                ax.scatter(df_subset[reference_region], df_subset[region], c=df_subset["diff"], cmap=cmap, alpha=0.7, edgecolors="none", norm=norm, rasterized=True)
                ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=0.8)
                
                ax.set_xlim(-0.005, 0.255)
                ax.set_ylim(-0.005, 0.255)
                ax.set_xlabel("")  # Remove x and y labels
                ax.set_ylabel("")  # Already set above
                
                # Set ticks for better visualization
                nticks = 6
                maxticks = 0.25
                ax.set_xticks(np.linspace(0, maxticks, nticks))
                ax.set_yticks(np.linspace(0, maxticks, nticks))

                # Add letter to each subplot
                ax.text(0.03, 0.97, letters[j*n_serotypes + i], transform=ax.transAxes, fontsize=36, fontweight="bold", va="top", ha="left")
            

        # Add the "VP1 Distance" label for the entire figure
        fig.text(0.5, 0.015, "VP1 Distance", ha="center", fontsize=36)

        # Add column headers for each serptype and row headers for each region
        for ax, region in zip(axes[:, 0], compare_regions):
            ax.set_ylabel(f"{region} Distance", fontsize=36)

        for ax, serotype in zip(axes[0], pdistance_dict.keys()):
            ax.set_title(f"{serotype}", fontsize=40, weight="bold")

        # Leave space for the colorbar
        plt.tight_layout(rect=[0, 0.02, 1, 0.98])

        # Add colorbar
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, location="top", fraction=0.02, pad=0.04)
        cbar.set_label("Absolute Difference in Hamming Distance", fontsize=40, labelpad=20)
        cbar.ax.tick_params(labelsize=36)

        # Save the plot
        output_file = f"{outdir}/grid_pdistances.{args.outformat}"
        plt.savefig(output_file)
        

        
    else:
        
        # For each serotype, plot the pairwise p-distances between the genomic regions
        for serotype, df in pdistance_dict.items():
            
            # For each compare_regions, plot the pairwise p-distances against the reference_region
            for region in compare_regions:
                
                # Filter the dataframe to only include the two regions
                df_subset = df[["seq1", "seq2", reference_region, region]].copy()
                
                # Remove rows where with missing values
                df_subset = df_subset.dropna()
                
                # Add a new column with the difference between the two regions
                df_subset["diff"] = abs(df_subset[region] - df_subset[reference_region])
                
                # Make scatter plot
                plt.figure(figsize=(8, 8))
                sns.scatterplot(data=df_subset, x=reference_region, y=region, alpha=0.7, linewidth=0, hue="diff", palette="viridis", hue_norm=norm)
                plt.title(f"{serotype}: {reference_region} vs. {region}")
                plt.xlabel(f"{reference_region} p-distance")
                plt.ylabel(f"{region} p-distance")
                # Plot diagonal line
                sns.lineplot(x=[0, 1], y=[0, 1], color="black", linestyle="--", linewidth=1)
                # Limit axes from 0 to highest value
                max_val = max(df_subset[reference_region].max(), df_subset[region].max())
                plt.xlim(-0.01, 1.1*max_val)
                plt.ylim(-0.01, 1.1*max_val)
                plt.tight_layout()
                
                # Save the plot
                output_file = f"{outdir}/{serotype}_{reference_region}_vs_{region}_pdistances.{args.outformat}"
                plt.savefig(output_file)
                
                    
                
            
                  

