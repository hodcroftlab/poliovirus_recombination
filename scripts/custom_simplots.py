### Custom script for plotting SimPlots (similarity plots)

# Input: alignment file (fasta)
# Parameters: windowsize, stepsize, sequence_id of the reference
# Output: similarity plot

# Use coloring scheme I used for the trees

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import argparse 
import os

argparser = argparse.ArgumentParser(description="Create custom SimPlots.")
argparser.add_argument("--alignment", help="Alignment file (fasta).")
argparser.add_argument("--windowsize", help="Window size for SimPlots.", type=int)
argparser.add_argument("--stepsize", help="Step size for SimPlots.", type=int)
argparser.add_argument("--metadata", help="Metadata file (CSV/TSV) containing serotype information.")
argparser.add_argument("--colors", help="CSV/TSV file containing colors for serotypes.")
argparser.add_argument("--exclude_neighbors", action="store_true", help="Exclude neighbors from the similarity plot.")
argparser.add_argument("--sabin", action="store_true", help="Replace PV consensus sequences with Sabin vaccine strain sequences.")
argparser.add_argument("--outdir", help="Output directory for SimPlots.")
argparser.add_argument("--outformat", help="Output format for the plots (e.g., svg, png).")
args = argparser.parse_args()

# Set font to Arial
plt.rcParams["font.family"] = "Arial"

# Function to split the alignment into windows of the given window size and step size
def split_alignment(alignment, window_size, step_size):
    windows = {} # Initialize a dictionary to store the windows
    sequence_length = len(alignment[0].seq)
    
    # Split the alignment into windows of the given window size, centered around the step size * i
    # At edges, the window size will be smaller
    for i in range(1, len(alignment[0].seq) // step_size + 1):
        center = i * step_size
        start = center - window_size // 2
        end = center + window_size // 2
        
        if start < 0:
            start = 0
        
        if end > sequence_length:
            end = sequence_length
            
        # Initialize alignment for each window
        window_alignment = []
            
        for record in alignment:
            # Slice the sequence and make sure it's a Seq object
            sub_seq = record.seq[start:end]
            
            # Create a new SeqRecord object with the sliced sequence
            new_record = record[:]  # Make a shallow copy of the record
            new_record.seq = sub_seq  # Assign the sliced sequence
            
            # Append the new record to the window alignment
            window_alignment.append(new_record)
            
        # Store the windowed alignment at the center position
        windows[center] = window_alignment
        
    return windows


# Function to calculate pairwise distances between a reference sequence and all other sequences in the alignment
def calculate_pairwise_distances(alignment, reference_id, current_step):

    # Get reference sequence from alignment
    reference_seq = [record for record in alignment if record.id == reference_id][0]
    reference_seq = np.array(list(reference_seq))
    
    # Remove the reference sequence from the alignment
    alignment = [record for record in alignment if record.id != reference_id]
    
    # Intialize results list
    results = []
    
    for record in alignment:
        query_seq = np.array(list(record.seq))
            
        # Mask all positions where either of the two sequences has a gap or an N
        valid_positions = (reference_seq != "-") & (reference_seq != "N") & (query_seq != "-") & (query_seq != "N")
        reference_valid = reference_seq[valid_positions]
        query_valid = query_seq[valid_positions]
        
        # Update the sequence length
        seq_len_valid = len(reference_valid)
        
        if seq_len_valid < args.windowsize/10:
            print(f"Skipping step {current_step} for {record.id} because it has fewer than 10% valid positions.")
            continue
        
        # Now calculate the number of differing positions
        nd = np.sum(query_valid != reference_valid)
        
        # Calculate the p-distance
        pdist = nd / seq_len_valid
        pdist = round(pdist, 4) # Round to 4 decimal places
        
        # Append the result as a tuple to the results list
        results.append((record.id, reference_id, pdist))
            
    # Convert the list of results to a dataframe
    results_df = pd.DataFrame(results, columns=["seq1", "seq2", "pdist"])
    
    # Add similarity column
    results_df["similarity"] = 1 - results_df["pdist"]
    
    return results_df

# Function to set colors for the sequences based on serotype information
def set_colors(results_df, metadata, colors):
    
    # Find consensus sequences as seq1 values that do not contain a dot (".")
    consensus_seqs = results_df["seq1"].unique()
    consensus_seqs = [seq for seq in consensus_seqs if "." not in seq]
    consensus_seqs = [seq for seq in consensus_seqs if seq not in ["EVC", "non-polio"]] # Remove "EVC" and "non-polio" from consensus sequences
    
    # Add consensus to metadata
    metadata = metadata[["Accession", "serotype"]]
    metadata = pd.concat([metadata, pd.DataFrame({"Accession": consensus_seqs, "serotype": consensus_seqs})])
    
    # Merge with metadata for serotype information
    results_df = results_df.merge(metadata[["Accession", "serotype"]], left_on="seq1", right_on="Accession")
    results_df.rename(columns={"serotype": "serotype1"}, inplace=True)
    
    # Merge again for seq2 (query sequence)
    results_df = results_df.merge(metadata[["Accession", "serotype"]], left_on="seq2", right_on="Accession")
    results_df.rename(columns={"serotype": "serotype2"}, inplace=True)
    
    # Merge with colors for color information
    results_df = results_df.merge(colors, left_on="serotype1", right_on="value")
    #results_df.drop(["class", "value", "Accession"], axis=1, inplace=True)
    
    return results_df

# Function to generate and save the SimPlots
def plot_simplot(results_df, outdir, outformat):
    
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
    elif "PV1" in sequences:
        pv_seqs = [seq for seq in sequences if "PV" in seq]
        other_seqs = [seq for seq in sequences if seq not in pv_seqs] # Get non-pv sequences and sort alphabetically
        other_seqs = sorted(other_seqs)
        sequences = pv_seqs + other_seqs
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
    
    ax.set_title(f"Query Sequence: {reference_seq}", fontsize=28)
    ax.set_xlabel("Position", fontsize=24)
    ax.set_ylabel("Similarity", fontsize=24)
    ax.legend(
        loc="center",                   # Center the legend box itself
        bbox_to_anchor=(1.18, 0.5),     # Position relative to the right side of the plot
        fontsize=12,
        ncol=2,
        borderaxespad=0,                # Space between the legend and the axes
        frameon=False                   # Remove the frame for a cleaner look
    )
    # Set fontsize of tick labels
    ax.tick_params(axis="both", which="major", labelsize=20)
    ax.set_ylim(0.4, 1.02)
    reference_seq = reference_seq.replace(" ", "_").replace("(", "").replace(")", "")
    plt.tight_layout()
    plt.savefig(f"{outdir}/simplot_{reference_seq}.{outformat}")
    ax.clear()
    plt.close(fig)
    
    
def main():
    
    # Read alignment
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    sequence_id = args.alignment.split("/")[-1].split(".fas")[0] # ID of the reference sequence
    
    # Make Sabin dict for later renaming
    sabin_dict = {"AY184219.1": "OPV1", "AY184220.1": "OPV2", "AY184221.1": "OPV3"}
    
    # Adjust data depending on exclude_neighbors and sabin options
    if args.sabin == True: # Remove PV consensus sequences from the data
        pv_consensus = ["PV1", "PV2", "PV3"]
        alignment = [record for record in alignment if record.id not in pv_consensus]
    else: # Remove Sabin sequences from the data
        alignment = [record for record in alignment if record.id not in sabin_dict.keys()]
            
    if args.exclude_neighbors == True: # Remove neighbor sequences from the data (identified by "."), except the reference sequence and the Sabin sequences
        neighbor_seqs = [record.id for record in alignment if "." in record.id]
        # Remove the reference sequence and Sabin sequences from the list
        neighbor_seqs = [seq_id for seq_id in neighbor_seqs if seq_id not in sabin_dict.keys()]
        neighbor_seqs.remove(sequence_id)
        alignment = [record for record in alignment if record.id not in neighbor_seqs]
    
    
    # Split the alignment into windows
    windows = split_alignment(alignment, args.windowsize, args.stepsize)
    
    # Initialize a dataframe with columns step, seq1, seq2, pdist to store the results
    results = pd.DataFrame(columns=["step", "seq1", "seq2", "pdist", "similarity"])
    
    # Calculate pairwise distances for each window
    for step, alignment in windows.items():
        window_results = calculate_pairwise_distances(alignment, sequence_id, current_step=step)
        window_results["step"] = step
        results = pd.concat([results, window_results])
        
    # Read metadata and colors
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
    else:
        raise ValueError("Please provide a CSV or TSV file for the metadata.")
    
    if args.colors.endswith(".csv"):
        colors = pd.read_csv(args.colors, names=["class", "value", "color"])
    elif args.colors.endswith(".tsv"):
        colors = pd.read_csv(args.colors, sep="\t", names=["class", "value", "color"])
    else:
        raise ValueError("Please provide a CSV or TSV file for the colors.")

    # Use abbreviated serotype in metadata and colors
    metadata["serotype"] = metadata["serotype"].str.replace("Enterovirus C", "EVC")
    metadata["serotype"] = metadata["serotype"].str.replace("Coxsackievirus A", "CVA")
    metadata["serotype"] = metadata["serotype"].str.replace("Poliovirus ", "PV")

    colors["value"] = colors["value"].str.replace("Enterovirus C", "EVC")
    colors["value"] = colors["value"].str.replace("Coxsackievirus A", "CVA")
    colors["value"] = colors["value"].str.replace("Poliovirus ", "PV")

    
    # Set colors for the sequences
    results = set_colors(results, metadata, colors)
    
    # Rename Sabin sequences in results
    results["seq1"] = results["seq1"].replace(sabin_dict)
    
    # Modify the output directory based on the exclude_neighbors and sabin options
    if args.exclude_neighbors == True:
        if args.sabin == True:
            outdir = f"{args.outdir}/without_neighbors/sabin"
        else:
            outdir = f"{args.outdir}/without_neighbors/no_sabin"
    else:
        if args.sabin == True:
            outdir = f"{args.outdir}/with_neighbors/sabin"
        else:
            outdir = f"{args.outdir}/with_neighbors/no_sabin"
            
    # Check if outdir exists, if not create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Plot the SimPlots
    plot_simplot(results, outdir, args.outformat)
    
if __name__ == "__main__":
    main()