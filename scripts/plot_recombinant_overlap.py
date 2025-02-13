### This script makes venn diagrams of the overlap between recombinant sequences identified via SimPlots and HMM

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import argparse

if __name__ == "__main__":
    
    argparser = argparse.ArgumentParser(description="Make venn diagrams of the overlap between recombinant sequences identified via SimPlots and HMM.")
    argparser.add_argument("--simplot_metadata", help="Metadata file containing recombination status inferred by SimPlots for each sequence.")
    argparser.add_argument("--hmm_results", help="Array containing the results of the HMM run.")
    argparser.add_argument("--outcsv", help="Output path for CSV of identified breakpoint locations by HMM.")
    argparser.add_argument("--outplot", help="Output path for venn diagram.")

    args = argparser.parse_args()
    
    # Read SimPlots metadata
    #simplot_metadata = pd.read_csv("../data/metadata/evc_full_genomes_metadata_recombination_status.csv")
    simplot_metadata = pd.read_csv(args.simplot_metadata)
    total_sequences = len(simplot_metadata)
    recombinant_sequences_simplot = simplot_metadata[simplot_metadata["recombination_status"] == "recombinant"]
    
    # Read HMM results
    #hmm_results = pd.read_csv("../data/recombination_analysis/hmm_recombination/hmm_results.csv")
    hmm_results = pd.read_csv(args.hmm_results)
    
    
    ### Extract recombination status and breakpoints from HMM results
    
    # Iterate over rows in the HMM dataframe
    # For each row, find the nucleotide positions at which the serotypes change
    breakpoints_df = pd.DataFrame(columns=["sequence", "breakpoint"])
    for idx, row in hmm_results.iterrows():
        inferred_serotypes = list(row[3:])
        breakpoints = []
        # Find indices where the serotype changes
        for i in range(len(inferred_serotypes)-1):
            if inferred_serotypes[i] != inferred_serotypes[i+1]:
                breakpoints.append(i)
                
        # Save the breakpoints for each sequence
        breakpoints_df = pd.concat([breakpoints_df, pd.DataFrame({"sequence": row["accession"], "breakpoint": breakpoints})])

    # Save the breakpoints to a CSV file
    #breakpoints_df.to_csv("../data/recombination_analysis/hmm_recombination/hmm_breakpoints.csv", index=False)
    breakpoints_df.to_csv(args.outcsv, index=False)
    
    
    ### Prepare data for venn diagram
    hmm_recombinants = set(breakpoints_df["sequence"].unique())
    simplot_recombinants = set(recombinant_sequences_simplot["accession"])
    
    # Find the overlap between the two sets
    overlap = hmm_recombinants.intersection(simplot_recombinants)
    
    # Find sequences that are unique to HMM
    unique_hmm = hmm_recombinants - simplot_recombinants
    # Check how many of these sequences were classified as unclear by SimPlot analysis
    unclear_sequences = simplot_metadata[simplot_metadata["recombination_status"] == "unclear"]
    unclear_sequences = set(unclear_sequences["accession"])
    unclear_hmm = unique_hmm.intersection(unclear_sequences)
    percentage = len(unclear_hmm) / len(unique_hmm) * 100
    print(f"Percentage of sequences identified as recombinant by HMM and unclear by SimPlot: {percentage:.0f}%")
    
    
    ### Create the venn diagram
    
    # Create the plot
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 30
    })
    fig, ax = plt.subplots(figsize=(6, 6))
    
    venn = venn2([hmm_recombinants, simplot_recombinants], set_labels = ("HMM", "SimPlot"), set_colors=("#005AB5", "#DC3220"))
    venn.get_label_by_id("A").set_color("#005AB5")
    venn.get_label_by_id("B").set_color("#DC3220")
    venn.get_label_by_id("10").set_fontweight("bold")
    venn.get_label_by_id("01").set_fontweight("bold")
    venn.get_label_by_id("11").set_fontweight("bold")
    
    # Save the plot
    plt.savefig(args.outplot, bbox_inches="tight")