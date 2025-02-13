### Custom script for plotting SimPlots (similarity plots)

# Input: alignment file (fasta)
# Parameters: windowsize, stepsize, sequence_id of the reference
# Output: similarity plot

# Use coloring scheme I used for the trees

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import argparse 
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("agg") # Use the Agg backend to save plots without displaying them

argparser = argparse.ArgumentParser(description="Create custom SimPlots.")
argparser.add_argument("--alignment", help="Fasta file of aligned sequences; make a SimPlot for each sequence in this alignment.")
argparser.add_argument("--neighbors", help="CSV/TSV containing neighbors to be plotted for each sequece.")
argparser.add_argument("--consensus", help="Fasta file of consensus sequences.")
argparser.add_argument("--windowsize", help="Window size for SimPlots.", type=int)
argparser.add_argument("--stepsize", help="Step size for SimPlots.", type=int)
argparser.add_argument("--metadata", help="Metadata file (CSV/TSV) containing serotype information.")
argparser.add_argument("--colors", help="CSV/TSV file containing colors for serotypes.")
argparser.add_argument("--outdir_plot", help="Output directory for SimPlots.")
argparser.add_argument("--outformat", help="Output format for the plots (e.g., png).")
argparser.add_argument("--outdir_csv", help="Output directory for csv files with similarity results for each query.")
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
    reference_seq = np.array(list(reference_seq.seq))
    
    # Remove the reference sequence from the alignment
    alignment = [record for record in alignment if record.id != reference_id]
    
    # Intialize results list
    results = []
    
    for record in alignment:
        query_seq = np.array(list(record.seq))
            
        # Mask all positions where either of the two sequences has a gap or an N
        valid_positions = (reference_seq != "-") & (reference_seq != "N") & (reference_seq != "n") & (query_seq != "-") & (query_seq != "N") & (query_seq != "n")
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
    # Set fontsize of tick labels
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_ylim(0.4, 1.02)
    reference_seq = reference_seq.replace(" ", "_").replace("(", "").replace(")", "")
    plt.tight_layout()
    plt.savefig(f"{outdir}/simplot_{reference_seq}.{outformat}")
    ax.clear()
    plt.close(fig)
    
    
    
def main():
    
    # Read alignment
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    
    # Convert all ambiguous nucleotides to N
    allowed_characters = ["A", "a", "T", "t", "G", "g", "C", "c", "N", "n", "-"]
    for record in alignment:
        # Convert all characters not in allowed_characters to N
        cleaned_seq = "".join([char if char in allowed_characters else "N" for char in record.seq])
        record.seq = Seq(cleaned_seq)
        
    # Get Sabin sequences
    # Dict for fetchig sequences and later renaming
    sabin_dict = {"AY184219.1": "OPV1", "AY184220.1": "OPV2", "AY184221.1": "OPV3"}
    sabin_records = [record for record in alignment if record.id in sabin_dict.keys()]
    
    # Remove from alignment
    alignment = [record for record in alignment if record.id not in sabin_dict.keys()]
    # This is the reason why input is 1940 sequences, but output is 1937 SimPlots
    
    # Get consensus sequences
    consensus = list(SeqIO.parse(args.consensus, "fasta"))
    consensus = [record for record in consensus if record.id not in ["EVC", "non-polio", "PV1", "PV2", "PV3"]] # Remove unwanted consensus sequences
    
    # Read neighbor data
    if args.neighbors.endswith(".csv"):
        neighbor_df = pd.read_csv(args.neighbors)
    elif args.neighbors.endswith(".tsv"):
        neighbor_df = pd.read_csv(args.neighbors, sep="\t")
    else:
        raise ValueError("Please provide a CSV or TSV file for the neighbors.")

    
    # Read metadata and colors
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
    else:
        raise ValueError("Please provide a CSV or TSV file for the metadata.")
    
    metadata = metadata[["Accession", "serotype"]] # select relevant columns

    if args.colors.endswith(".csv"):
        colors = pd.read_csv(args.colors, names=["class", "serotype", "color"])
    elif args.colors.endswith(".tsv"):
        colors = pd.read_csv(args.colors, sep="\t", names=["class", "serotype", "color"])
    else:
        raise ValueError("Please provide a CSV or TSV file for the colors.")
    
    colors = colors[["serotype", "color"]]
    
    # Use abbreviated serotypes in metadata and colors
    metadata["serotype"] = metadata["serotype"].str.replace("Enterovirus C", "EVC")
    metadata["serotype"] = metadata["serotype"].str.replace("Coxsackievirus A", "CVA")
    metadata["serotype"] = metadata["serotype"].str.replace("Poliovirus ", "PV")
    metadata["serotype"] = metadata["serotype"].str.replace("Poliovirus", "PV") # For unidentified PV

    colors["serotype"] = colors["serotype"].str.replace("Enterovirus C", "EVC")
    colors["serotype"] = colors["serotype"].str.replace("Coxsackievirus A", "CVA")
    colors["serotype"] = colors["serotype"].str.replace("Poliovirus ", "PV")
    colors["serotype"] = colors["serotype"].str.replace("Poliovirus", "PV") # For unidentified PV


    # Add consensus sequences to metadata
    consensus_seqs = [record.id for record in consensus]
    metadata = pd.concat([metadata, pd.DataFrame({"Accession": consensus_seqs, "serotype": consensus_seqs})])
    
    # Merge metadata and colors to a single dataframe
    metadata_colors = metadata.merge(colors, on="serotype")
    
    
    # Iterate over all sequences in the alignment
    for record in alignment:
        sequence_id = record.id
        print(f"Creating SimPlot for {sequence_id}")
    
        # Get neigbors for the query sequence
        neighbor_list = list(set(neighbor_df[neighbor_df["seq"] == sequence_id]["neighbor_accn"].tolist())) # set() to remove duplicates
        
        # Fetch neighbor sequences from alignment
        neighbor_records = [record for record in alignment if record.id in neighbor_list]
        
        # Build alignment for the SimPlot
        simplot_alignment = [record for record in alignment if record.id == sequence_id]   # Add query sequence
        simplot_alignment += neighbor_records   # Add neighbor sequences
        simplot_alignment += consensus   # Add consensus sequences
        simplot_alignment += sabin_records   # Add Sabin sequences
        
        
        # Split the alignment into windows
        windows = split_alignment(simplot_alignment, args.windowsize, args.stepsize)
    
        # Initialize a dataframe with columns step, seq1, seq2, pdist to store the results
        results = pd.DataFrame(columns=["step", "seq1", "seq2", "pdist", "similarity"])
        
        # Calculate pairwise distances for each window
        for step, aln in windows.items():
            window_results = calculate_pairwise_distances(alignment=aln, reference_id=sequence_id, current_step=step)
            window_results["step"] = step
            results = pd.concat([results, window_results])
    
        # Set colors for the sequences by merging with metadata_colors
        # Check if query sequence is in metadata_colors
        # If not, skip the sequence	
        if sequence_id not in metadata_colors["Accession"].values:
            print(f"Skipping {sequence_id} due to missing metadata (usually caused by missing / too gappy VP1 sequence)")
            continue
        
        results = results.merge(metadata_colors, left_on="seq1", right_on="Accession", how="left")
        results = results.rename(columns={"serotype": "serotype1"})
        results = results.drop(columns=["Accession"])
        
        results = results.merge(metadata, left_on="seq2", right_on="Accession", how="left")
        results = results.rename(columns={"serotype": "serotype2"})
        results = results.drop(columns=["Accession"])
        
        
        # Rename Sabin sequences in results
        results["seq1"] = results["seq1"].replace(sabin_dict)
            
        # Check if output directories exists, if not create them
        if not os.path.exists(args.outdir_plot):
            os.makedirs(args.outdir_plot)
        if not os.path.exists(args.outdir_csv):
            os.makedirs(args.outdir_csv)
            
        # Save the results as a CSV file
        results.to_csv(f"{args.outdir_csv}/{sequence_id}_similarity_results.csv", index=False)
    
        # Plot the SimPlots
        plot_simplot(results, args.outdir_plot, args.outformat)
    
if __name__ == "__main__":
    main()