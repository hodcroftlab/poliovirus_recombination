### Re-plot SimPlots integrating the inferred breakpoints

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

argparser = argparse.ArgumentParser(description="Re-plot SimPlots integrating the inferred breakpoints.")
argparser.add_argument("--similarity_data", nargs="+", help="CSV file containing pairwise sequence similarity results.")
argparser.add_argument("--breakpoints_csv", help="CSV file containing inferred breakpoints.")
argparser.add_argument("--colors", help="CSV or TSV file containing colors for each serotype.")
argparser.add_argument("--outdir_plot", help="Output directory for the re-plotted SimPlots.")
argparser.add_argument("--outformat", help="Output format for the re-plotted SimPlots (png/svg/pdf).")
args = argparser.parse_args()

# Set global font style and size#
plt.rcParams.update({
    "font.family": "Arial",   
    "font.size": 20            
})

# Create output directory if it does not exist
if not os.path.exists(args.outdir_plot):
    os.makedirs(args.outdir_plot)

# Function to generate and save the SimPlots
def plot_simplot(results_df, breakpoints, outdir, outformat):
    
    # Add serotype information to the non-consensus sequences
    results_df["seq2"] = results_df["seq2"] + " (" + results_df["serotype2"] + ")" # Query sequence
    for seq in results_df["seq1"].unique():
        if "." in seq:
            results_df.loc[results_df["seq1"] == seq, "seq1"] = seq + " (" + results_df[results_df["seq1"] == seq]["serotype1"].values[0] + ")"
    
    # Make the similarity plot
    fig, ax = plt.subplots(figsize=(18, 5))
    
    # Get unique sequences
    sequences = results_df["seq1"].unique()
    
    # Set order of sequences for plotting: OPV1, OPV2, OPV3 / PV1, PV2, PV3; other serotypes
    if "OPV1" in sequences:
        sabin_seqs = [seq for seq in sequences if "OPV" in seq]
        other_seqs = [seq for seq in sequences if seq not in sabin_seqs] # Get non-sabin sequences and sort alphabetically
        other_seqs = sorted(other_seqs)
        sequences = sabin_seqs + other_seqs
    else:
        sequences = sorted(sequences)

    # Get reference sequence
    reference_seq = results_df["seq2"].unique()[0]
    
    # Plot each sequence
    for seq in sequences:
        # Plot neighbor sequences (identified by "." in the sequence name) as dashed lines
        # Make sure that if a step does not exist, the gap is correctly plotted
        if "." in seq:
            seq_results = results_df[results_df["seq1"] == seq]
            ax.plot(seq_results["step"], seq_results["similarity"], label=seq, color=seq_results["color"].values[0], linestyle="--", linewidth=1)
        else:
            seq_results = results_df[results_df["seq1"] == seq]
            ax.plot(seq_results["step"], seq_results["similarity"], label=seq, color=seq_results["color"].values[0])
    
    ax.set_title(f"Query Sequence: {reference_seq}", fontsize=20)
    ax.set_xlabel("Position", fontsize=20)
    ax.set_ylabel("Similarity", fontsize=20)
    ax.legend(
        loc="center",                   # Center the legend box itself
        bbox_to_anchor=(1.18, 0.5),     # Position relative to the right side of the plot
        fontsize=12,
        ncol=2,
        borderaxespad=0,                # Space between the legend and the axes
        frameon=False                   # Remove the frame for a cleaner look
    )
    # Add vertical black line at each breakpoint
    for bp in breakpoints:
        ax.axvline(x=bp, color="black", linestyle="--", linewidth=1)
    # Set fontsize of tick labels
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_ylim(0.4, 1.02)
    reference_seq = reference_seq.replace(" ", "_").replace("(", "").replace(")", "")
    plt.tight_layout()
    plt.savefig(f"{outdir}/simplot_with_bp_{reference_seq}.{outformat}")
    ax.clear()
    plt.close(fig)
    
# Read breakpoints
breakpoints = pd.read_csv(args.breakpoints_csv)
    
# Read data from the similarity results files
for csv in args.similarity_data:
    
    sequence = csv.split("/")[-1].split("_")[0]
    data = pd.read_csv(csv)
    
    print("Re-plotting SimPlots for", sequence, "...")
    
    # Get list of breakpoints from the breakpoints dataframe
    specific_breakpoints = breakpoints[breakpoints["sequence"] == sequence]
    specific_breakpoints = specific_breakpoints["breakpoint"].tolist()
    
    print("Breakpoints for", sequence, ":", specific_breakpoints)
    
    # Make the similarity plot
    plot_simplot(data, specific_breakpoints, args.outdir_plot, args.outformat)