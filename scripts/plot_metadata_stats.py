### This scripts generates plots showing some statistics of the sequencing data (e.g., number of sequences per year and continent, )

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
import os

if __name__ == "__main__":
    
    # Parse arguments
    argparser = argparse.ArgumentParser(description="Plot number of sequences per serotype.")
    argparser.add_argument("--metadata", help="Metadata file (CSV/TSV).")
    argparser.add_argument("--colors", help="CSV/TSV file containing mapping between serotypes and colors.")
    argparser.add_argument("--outdir", help="Output directory for the plot.")
    
    args = argparser.parse_args()
    
    # Check if output directory exists, if not, create it
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Set plotting parameters
    plt.rcParams.update({"font.size": 14})
    plt.rcParams.update({"axes.labelsize": 16})
    plt.rcParams.update({"xtick.labelsize": 14})
    plt.rcParams.update({"ytick.labelsize": 14})
    plt.rcParams.update({"legend.fontsize": 14})
    plt.rcParams["font.family"] = "Arial"
    
    # Read metadata
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
    else:
        raise ValueError("Metadata file must be in CSV or TSV format.")


    ## Plot pie chart and bar chart of the number of sequences for each serotype
    
    # Count sequences per serotype
    df_serotypes = metadata.groupby("serotype_short").size().reset_index(name="count")
    df_serotypes.rename(columns={"serotype_short": "serotype"}, inplace=True)
    
    # Group serotypes with less than 10 sequences as "Other"
    df_serotypes.loc[df_serotypes["count"] < 10, "serotype"] = "Other"
    df_serotypes = df_serotypes.groupby("serotype").sum().reset_index()
    
    # Sort by total number of sequences, but keep "Other" at the end
    df_serotypes = df_serotypes.sort_values(by="count", ascending=False)
    df_serotypes = pd.concat([df_serotypes[df_serotypes["serotype"] != "Other"], df_serotypes[df_serotypes["serotype"] == "Other"]]) # Append "Other" at the end using pd.concat
    
    # Get color mapping for serotypes
    if args.colors.endswith(".tsv"):
        colors = pd.read_csv(args.colors, sep="\t", names=["class", "serotype", "color"])
    elif args.colors.endswith(".csv"):
        colors = pd.read_csv(args.colors, names=["class", "serotype", "color"])
    else:
        raise ValueError("Color file must be in CSV or TSV format.")
    
    colors["serotype"] = colors["serotype"].str.replace("Enterovirus C", "EVC")
    colors["serotype"] = colors["serotype"].str.replace("Coxsackievirus A", "CVA")
    colors["serotype"] = colors["serotype"].str.replace("Poliovirus ", "PV")
    
    # Add color for "Other" serotype
    colors = pd.concat([colors, pd.DataFrame({"class": ["serotype"], "serotype": ["Other"], "color": ["#808080"]})], ignore_index=True)
    
    # Merge colors with serotype counts
    df_serotypes = df_serotypes.merge(colors, on="serotype", how="left")
    
    # Make pie chart showing counts of full genomes per serotype (not percentage)
    # Custom autopct function to reduce the font size of counts
    def custom_autopct(pct):
        total = values.sum()
        val = int(pct * total / 100)
        return '{:.0f}'.format(val)
    
    values = df_serotypes["count"]
    plt.figure(figsize=(6, 6))
    wedges, texts, autotexts = plt.pie(
        df_serotypes["count"], 
        labels=df_serotypes["serotype"], 
        startangle=0, 
        autopct=custom_autopct,  # Use custom function for counts
        colors=df_serotypes["color"],  # Use custom colors
        pctdistance=0.8,  # Adjust distance of counts
        textprops={"fontsize": 11},  # Set font size for labels
    )

    # Change the font size of the autotexts (counts)
    for autotext in autotexts:
        autotext.set_fontsize(8)  # Set the font size for counts
    
    plt.title("Nr. of Full Genome Sequences")
    
    # Save plot as .png
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/piechart_serotypes_full_genomes.png", dpi=300)
    
    # Now a bar chart showing the same data
    df_serotypes.plot(x="serotype", y="count", kind="bar", width=0.9, color=df_serotypes["color"], legend=False, figsize=(8, 6))
    plt.xlabel("Serotype")
    plt.ylabel("Number of Sequences")
    plt.title("Number of Full Genome Sequences per Serotype")
    
    # Extend y-axis by 10% to make space for text labels
    plt.ylim(0, df_serotypes["count"].max() * 1.1)
    
    # Add total number of sequences as text on top of the bars
    for i in range(len(df_serotypes)):
        plt.text(x=i, y=df_serotypes.iloc[i]["count"]+15, s=int(df_serotypes.iloc[i]["count"]), ha="center")
        
    # Save plot as .png
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/barchart_serotypes_full_genomes.png", dpi=300)
    
    
    # Now a single, horizontal stacked bar of the same data
    # Labels with serotype and count within the bars
    
    # Example serotype data (replace with your df_serotypes data)
    serotypes = df_serotypes["serotype"]
    values = df_serotypes["count"]

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 4))

    # Initialize the starting position of the left edge of the bars
    left = 0
    others = 0

    # Loop through each serotype and value to create a stacked bar chart
    for i, (serotype, value) in enumerate(zip(serotypes, values)):
        # Plot each section of the stacked bar
        color = df_serotypes[df_serotypes["serotype"] == serotype]["color"].values[0]
        ax.barh(0, value, left=left, color=color, edgecolor="white")
        
        # If the value is small, sum the remaining values and plot them as "Other" below the stacked bar
        if value < 200 and "PV" not in serotype:
            others += value
        else:
            # Add text labels in the middle of each section
            ax.text(left + value / 2, 0, f"{serotype}\n{value}", va="center", ha="center", fontsize=24, color="black")
        
        # Update the left position for the next bar
        left += value

    # Remove axes (only show the stacked bar)
    ax.set_axis_off()
    
    # Plot total number of sequences above the stacked bar
    total = values.sum()
    ax.text(0, 0.5, f"Total Sequences: {total}", va="center", ha="left", fontsize=28, color="black")
    #ax.text(0, -0.5, f"Total Sequences: {total}", va="center", ha="left", fontsize=28, color="black")
    
    # Add "others" label below the stacked bar aligned to the right
    ax.text(left, -0.5, f"Non-Polio EVC: {others}", va="center", ha="right", fontsize=24, color="black")

    # Save plot as .png
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/stackedbar_serotypes_full_genomes.svg")

    
    
    
    