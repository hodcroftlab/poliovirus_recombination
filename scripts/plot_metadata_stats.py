### This scripts generates plots showing some statistics of the sequencing data (e.g., number of sequences per year and continent, )

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

if __name__ == "__main__":
    
    # Parse arguments
    argparser = argparse.ArgumentParser(description="Plot statistics of the sequencing data.")
    argparser.add_argument("--metadata_polio", help="Metadata file for polio sequences in csv format.")
    argparser.add_argument("--metadata_npevc", help="Metadata file for non-polio EVC sequences in csv format.")
    argparser.add_argument("--output_dir", help="Output directory for the plots.")
    
    args = argparser.parse_args()
    
    # Set plotting parameters
    plt.rcParams.update({"font.size": 12})
    plt.rcParams.update({"axes.labelsize": 14})
    plt.rcParams.update({"xtick.labelsize": 12})
    plt.rcParams.update({"ytick.labelsize": 12})
    plt.rcParams.update({"legend.fontsize": 12})
    plt.rcParams["font.family"] = "Arial"
    
    # Read metadata
    metadata_polio = pd.read_csv(args.metadata_polio)
    metadata_npevc = pd.read_csv(args.metadata_npevc)
    
    metadata_polio = pd.read_csv("../data/metadata/polio_metadata.csv")
    metadata_npevc = pd.read_csv("../data/metadata/npevc_metadata.csv")
    
    ## Merge metadata to plot bar chart of the number of sequences per serotype, stacked by full genome and partial sequences
    metadata = pd.concat([metadata_polio, metadata_npevc])
    
    # Edit Nuc_Completeness column by assigning "full" to all sequences longer than 6000 bp and "partial" to all others
    metadata["Nuc_Completeness"] = ["full" if x >= 6000 else "partial" for x in metadata["Length"]]
    
    df_serotypes = metadata.groupby(["serotype", "Nuc_Completeness"]).size().reset_index(name="count") # Count sequences per serotype and completeness
    
    # Replace "Enterovirus " by "EV-", "Coxsackievirus " by "CV-", and Poliovirus " by "PV-" in serotype names
    df_serotypes["serotype"] = df_serotypes["serotype"].str.replace("Enterovirus ", "EV-", regex=False)
    df_serotypes["serotype"] = df_serotypes["serotype"].str.replace("Coxsackievirus ", "CV-", regex=False)
    df_serotypes["serotype"] = df_serotypes["serotype"].str.replace("Poliovirus ", "PV-", regex=False)
    df_serotypes["serotype"] = df_serotypes["serotype"].str.replace("Poliovirus", "PV", regex=False)
    df_serotypes["serotype"] = df_serotypes["serotype"].str.replace("Pv", "PV", regex=False)
    
    # Pivot table to have "full" and "partial" as columns
    df_serotypes = df_serotypes.pivot(index="serotype", columns="Nuc_Completeness", values="count").reset_index().fillna(0)
    
    # Make stacked barplot ordered by total number of sequences
    df_serotypes["total"] = df_serotypes["full"] + df_serotypes["partial"]
    
    # Group serotypes with less than 10 total sequences as "Other"
    df_serotypes.loc[df_serotypes["total"] < 10, "serotype"] = "Other"
    df_serotypes = df_serotypes.groupby("serotype").sum().reset_index()
    
    # Sort by total number of sequences, but keep "Other" at the end
    df_serotypes = df_serotypes.sort_values(by="total", ascending=False)
    df_serotypes = pd.concat([df_serotypes[df_serotypes["serotype"] != "Other"], df_serotypes[df_serotypes["serotype"] == "Other"]]) # Append "Other" at the end using pd.concat
    
    df_serotypes = df_serotypes.drop(columns="total") # Drop "total" column
    
    # Make stacked barplot of full and partial sequences per serotype
    df_serotypes.plot(x="serotype", kind="bar", stacked=True, figsize=(12, 6), width=0.9, color=["#006BA4", "#FF800E"]) 
    plt.xlabel("Serotype")
    plt.ylabel("Number of Sequences")
    plt.title("Number of Full and Partial Sequences per Serotype")
    plt.legend(title="Genome Status", labels=["Full (>= 6000 bp)", "Partial (< 6000 bp)"], loc="upper right") # Adjust legend
    
    # Add total number of sequences as text on top of the bars
    for i in range(len(df_serotypes)):
        plt.text(x=i, y=df_serotypes.iloc[i]["full"] + df_serotypes.iloc[i]["partial"]+15, s=int(df_serotypes.iloc[i]["full"] + df_serotypes.iloc[i]["partial"]), ha="center")
    
    # Save plot as .png
    plt.tight_layout()
    plt.savefig(f"{args.output_dir}/sequences_per_serotype.png", dpi=300)
    
    
    ## Plot number of sequences per year as stacked histogram colored by continent, in separate subfigures for polio and non-polio EV-C
    
    # Add release year column
    metadata_polio["release_year"] = pd.to_datetime(metadata_polio["Release_Date"]).dt.year
    metadata_npevc["release_year"] = pd.to_datetime(metadata_npevc["Release_Date"]).dt.year
    
    # Convert NA values to "Unknown" for continent
    metadata_polio["continent"] = metadata_polio["continent"].fillna("Unknown")
    metadata_npevc["continent"] = metadata_npevc["continent"].fillna("Unknown")
    
    # Keep only sequences starting in 2000
    metadata_polio = metadata_polio[metadata_polio["release_year"] >= 2000]
    metadata_npevc = metadata_npevc[metadata_npevc["release_year"] >= 2000]

    # Plot number of sequences per year as stacked barplot colored by continent
    fig, axs = plt.subplots(2, 1, figsize=(12, 12))
    fig.suptitle("Number of Sequences by Release Year and Continent", fontsize=16)
    
    df_polio = metadata_polio.groupby(["release_year", "continent"]).size().reset_index(name="count")
    df_npevc = metadata_npevc.groupby(["release_year", "continent"]).size().reset_index(name="count")
    
    # Pivot tables to have continents as columns
    df_polio = df_polio.pivot(index="release_year", columns="continent", values="count").reset_index().fillna(0)
    df_npevc = df_npevc.pivot(index="release_year", columns="continent", values="count").reset_index().fillna(0)

    # Make stacked barplots    
    df_polio.plot(x="release_year", kind="bar", stacked=True, ax=axs[0], color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999"], width=0.9)
    df_npevc.plot(x="release_year", kind="bar", stacked=True, ax=axs[1], color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999"], width=0.9)
    
    # Add labels, titles, and legends
    # Remove x-axis label for the first subplot
    axs[0].set_xlabel("")
    axs[0].set_ylabel("Number of Sequences")
    axs[0].set_title("Poliovirus", x=0.01, y=0.92, horizontalalignment="left")
    axs[0].get_legend().remove() # Remove legend for the first subplot
    axs[0].set_ylim(0, 1200) # Manually set y-axis limits to make plots comparable
    
    axs[1].set_xlabel("Release Year")
    axs[1].set_ylabel("Number of Sequences")
    axs[1].set_title("Non-Polio Enterovirus C", x=0.01, y=0.92, horizontalalignment="left")
    axs[1].legend(title="Continent", loc="center left")
    axs[1].set_ylim(0, 1200) # Manually set y-axis limits to make plots comparable
    
    # Save plot as .png
    plt.tight_layout()
    plt.savefig(f"{args.output_dir}/sequences_per_year_continent.png", dpi=300)
    
    
    
    
    
    
    
    
    