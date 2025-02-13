### Script for making SimPlots as subplots for the thesis 

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import argparse 
import os
import matplotlib.pyplot as plt

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
        
        if seq_len_valid < 200/10:
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
    results = pd.DataFrame(results, columns=["seq1", "seq2", "pdist"])
    
    # Add similarity column
    results["similarity"] = 1 - results["pdist"]
    
    return results
    
    
    
def main():
    
    # Read alignment
    alignment = list(SeqIO.parse("../data/sequences/alignments/full_genome/assembled/fullgenome_pre_masking_assembled.fasta", "fasta"))
    
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
    
    # Get consensus sequences
    consensus = list(SeqIO.parse("../data/sequences/consensus_sequences/all_consensus_sequences.fasta", "fasta"))
    consensus = [record for record in consensus if record.id not in ["EVC", "non-polio", "PV1", "PV2", "PV3"]] # Remove unwanted consensus sequences
    
    # Read neighbor data
    neighbor_df = pd.read_csv("../data/recombination_analysis/custom_simplots/neighbors_filtered_v2_all.csv")

    # Read metadata and colors
    metadata = pd.read_csv("../data/metadata/evc_full_genomes_metadata_rivm.csv")
    metadata = metadata[["Accession", "serotype"]] # Select relevant columns

    # Read colors
    colors = pd.read_csv("../data/config/colors.tsv", sep="\t", names=["class", "serotype", "color"])
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
    
    # Subset the alignment to sequences of interest
    sequences_of_interest = {
        "OPV-OPV": ["FJ460226.1", "FJ859184.1", "EU566935.1", "FJ859058.1", "AY948201.1"],
        "OPV-unknown": ["AF405688.1", "DQ890385.1", "HQ738299.1", "JX275031.2", "HF913427.1"],
        "OPV-nonpolio": ["AB180070.1", "AB180071.1", "AB180072.1", "AB180073.1", "MN149910.1"],
        "nonpolio-nonpolio": ["EU840733.1", "OK570217.1", "OK570238.1", "PP548242.1", "MZ171089.1"],
        "nonrecombinant": ["EU794954.1", "AY278549.1", "KJ170630.1", "AB769165.1", "KR815824.1"],
        "unclear": ["KJ019832.1", "OP137307.1", "OQ829421.1", "AB205396.1", "PP756352.1"]
    }
    
    # Make subalignments for each group of sequences
    subalignments = {}
    for group, seqs in sequences_of_interest.items():
        subalignments[group] = [record for seq_id in seqs for record in alignment if record.id == seq_id] # Preserve order of sequences
    
    
    # Iterate over subalignments
    for group, subalignment in subalignments.items():
        # Make one shared figure for each group, with a subplot for each sequence (2 columns, 3 rows, legend in 6th subplot)
        fig, axes = plt.subplots(3, 2, figsize=(18, 11))
        
        # Add a legend subplot instead of the last subplot
        fig.delaxes(axes[2, 1])
        legend_ax = fig.add_subplot(3, 2, 6)
        legend_ax.axis("off")

        # Add legend to the legend subplot
        for serotype, color in zip(colors["serotype"], colors["color"]):
            if serotype in ["Unidentified EVC", "Unidentified PV"]:
                continue
            elif "PV" in serotype:
                serotype = serotype.replace("PV", "OPV")
            legend_ax.plot([], [], label=serotype, color=color)

        # Add the legend with multiple columns
        legend = legend_ax.legend(
            title="Serotype", 
            loc="center", 
            fontsize=16, 
            title_fontsize=20,
            ncol=4  # Split legend into 4 columns
        )
        
        # Make the title bold
        legend.get_title().set_fontweight("bold")
        
        # Initialize letters for subplots
        letters = list("ABCDEF")
        
        for i, record in enumerate(subalignment):
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
            windows = split_alignment(simplot_alignment, 200, 50)
        
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
            
            # Add serotype information to the non-consensus sequences
            results["seq2"] = results["seq2"] + " (" + results["serotype2"] + ")" # Query sequence
            for seq in results["seq1"].unique():
                if "." in seq:
                    results.loc[results["seq1"] == seq, "seq1"] = seq + " (" + results[results["seq1"] == seq]["serotype1"].values[0] + ")"
            
            # Make the similarity plot
            ax = axes[i//2, i%2]
            
            # Get unique sequences
            sequences = results["seq1"].unique()
            
            # Set order of sequences for plotting: OPV1, OPV2, OPV3 / PV1, PV2, PV3; other serotypes
            if "OPV1" in sequences:
                sabin_seqs = [seq for seq in sequences if "OPV" in seq]
                other_seqs = [seq for seq in sequences if seq not in sabin_seqs] # Get non-sabin sequences and sort alphabetically
                other_seqs = sorted(other_seqs)
                sequences = sabin_seqs + other_seqs
            else:
                sequences = sorted(sequences)

            # Get reference sequence
            reference_seq = results["seq2"].unique()[0]
            
            # Store legend handles for neighbor sequences (to be plotted in extra legend in each subplot)
            neighbor_handles = []
            neighbor_labels = []
            
            # Plot each sequence
            for seq in sequences:
                
                seq_results = results[results["seq1"] == seq]
                
                # Plot neighbor sequences (identified by "." in the sequence name) as dashed lines
                # Make sure that if a step does not exist, the gap is correctly plotted
                if "." in seq:
                    line, = ax.plot(seq_results["step"], seq_results["similarity"], 
                                    label=seq, color=seq_results["color"].values[0], 
                                    linestyle="--", linewidth=1)
                    # Store the line and label for the neighbor sequences
                    neighbor_handles.append(line)
                    neighbor_labels.append(seq)
                else:
                    ax.plot(seq_results["step"], seq_results["similarity"], label=seq, color=seq_results["color"].values[0])
            
            # Add a separate legend in the bottom right for neighbor sequences (if any exist)
            if neighbor_handles:
                ax.legend(neighbor_handles, neighbor_labels, loc="lower right", 
                        fontsize=12, title_fontsize=14)
                    
            # Subplot letter
            subplot_letter = letters[i]
            ax.text(0.01, 0.98, subplot_letter, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top")
            
            ax.set_title(f"Query Sequence: {reference_seq}", fontsize=20)
            ax.set_xlabel("Position", fontsize=20)
            ax.set_ylabel("Similarity", fontsize=20)
            
            # Set fontsize of tick labels
            ax.tick_params(axis="both", which="major", labelsize=16)
            ax.set_ylim(0.38, 1.02)
            reference_seq = reference_seq.replace(" ", "_").replace("(", "").replace(")", "")
            
            # # Optional: add breakpoints
            # breakpoint_data = pd.read_csv("../data/recombination_analysis/custom_simplots/inferred_breakpoints.csv")
            # breakpoints = breakpoint_data[breakpoint_data["sequence"] == sequence_id]["breakpoint"].tolist()
            # # Add vertical black line at each breakpoint
            # for bp in breakpoints:
            #     ax.axvline(x=bp, color="black", linestyle="--", linewidth=1)
            
        plt.tight_layout(h_pad=2.0)
        plt.savefig(f"../plots/simplots/for_thesis/simplot_{group}.png")
                        
                
if __name__ == "__main__":
    main()