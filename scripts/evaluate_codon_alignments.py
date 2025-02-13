### This script evaluates the quality of a multiple sequence alignment of coding regions by counting the number of stop codons and incomplete codons in each sequence

from Bio import SeqIO
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description="Evaluate the quality of a multiple sequence alignment of coding regions by counting the number of stop codons and incomplete codons in each sequence.")
    argparser.add_argument("--alignments", help="Input alignment files to be evaluated.", nargs="+")
    argparser.add_argument("--output", help="Output file (csv) containing the total number of stop codons and incomplete codons in the alignment.")
    argparser.add_argument("--plot", help="Output file (png) containing a bar plot of the number of stop codons and incomplete codons in each alignment.")

    args = argparser.parse_args()

    alignment_files = args.alignments
    output_file = args.output

    # For each alignment file, for each sequence in the alignment, count the number of stop codons and incomplete codons
    # Also calculate average number of gaps and sum-of-pairs (SP) score
    # Output: dataframe with three columns: alignment file, number of sequences, length of the sequences, number of stop codons, number of incomplete codons
    
    results_df = pd.DataFrame(columns=["alignment", "nr_sequences", "sequence_length", "sp_score", "avg_gaps", "stop_codons", "incomplete_codons_ends", "incomplete_codons_internal"])
    
    def calculate_sp_score(alignment, match_score=1, mismatch_score=0, gap_penalty=-1):
        # Calculate the sum-of-pairs (SP) score for the alignment
        # The SP score is the sum of the scores for all pairs of sequences in the alignment
        # The score for a pair of sequences is determined by the number of matching nucleotides minus the number of mismatches minus the number of gaps
        
        # Convert alignment into a 2D NumPy array (rows = sequences, columns = alignment positions)
        sequences = np.array([list(seq.seq) for seq in alignment])  # Convert to character array
        n_sequences, n_columns = sequences.shape

        # Initialize SP score
        sp_score = 0

        # Iterate over each column
        for col_idx in range(n_columns):
            column = sequences[:, col_idx]  # Extract a single column
            unique_chars, counts = np.unique(column, return_counts=True)  # Count occurrences of each character

            # Precompute pairwise scores for this column
            for i, char1 in enumerate(unique_chars):
                for j, char2 in enumerate(unique_chars[i:], start=i):  # Avoid redundant pairs
                    count1, count2 = counts[i], counts[j]
                    pair_count = count1 * count2 if char1 != char2 else count1 * (count1 - 1) // 2

                    if char1 == char2:
                        if char1 == '-':  # Gap-gap match
                            sp_score += gap_penalty * pair_count
                        else:  # Nucleotide match
                            sp_score += match_score * pair_count
                    else:
                        if char1 == '-' or char2 == '-':  # Gap mismatch
                            sp_score += gap_penalty * pair_count
                        else:  # Nucleotide mismatch
                            sp_score += mismatch_score * pair_count

        return sp_score
    
    for alignment_file in alignment_files:
        alignment = list(SeqIO.parse(alignment_file, "fasta"))
        
        nr_sequences = len(alignment)
        sequence_length = len(alignment[0].seq)
        stop_codons_total = 0
        incomplete_codons_ends_total = 0
        incomplete_codons_internal_total = 0
        gaps_alignment = 0
        
        # Calculate the sum-of-pairs (SP) score for the alignment
        print(f"\nCalculating SP score for alignment {alignment_file}...")
        sp_score = calculate_sp_score(alignment)
          
        for record in alignment:
            # Check if all sequences within the alignment have the same length, otherwise the alignment is not valid
            if len(record.seq) != sequence_length:
                print(f"Sequence {record.id} in alignment {alignment_file} has a different length than the first sequence in the alignment (counting gaps). Make sure the sequences are aligned.")
                break
            
            # Alignments are nucleotide alignments, not AA alignments
            # Stop codons are TGA, TAA, TAG
            # Incomplete codons are sequences that contain 1 or 2 gap characters (e.g., "AT-", "A-G", "--T", "-A-""), but not 3 ("---", that would just be a gap in the sequence)
            # Differentiate between incomplete codons at the start and end of the sequence and in the middle of the sequence
            
            # Walk through the sequence in steps of 3
            stop_codons_seq = 0
            incomplete_codons_ends = 0
            incomplete_codons_internal = 0
            sequence = record.seq
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3]
                if codon in ["TGA", "TAA", "TAG", "tga", "taa", "tag"]:
                    stop_codons_seq += 1
                    print(f"Stop codon {codon} found starting at nucleotide position {i} in sequence {record.id} in alignment {alignment_file}.")
                elif "-" in codon and codon.count("-") < 3:
                    # Check whether the incomplete codon is at the start or end of the sequence, i.e., only gaps precede or follow the codon, respectively
                    if sequence[:i].count("-") == i or sequence[i+3:].count("-") == len(sequence) - i - 3:
                        incomplete_codons_ends += 1
                    else:
                        incomplete_codons_internal += 1
                        print(f"Internal incomplete codon {codon} found starting at nucleotide position {i} in sequence {record.id} in alignment {alignment_file}.")
            # Add the counts for the sequence to the alignment counts
            stop_codons_total += stop_codons_seq
            incomplete_codons_ends_total += incomplete_codons_ends
            incomplete_codons_internal_total += incomplete_codons_internal
            
            # Count the number of gaps in the sequence and add it to the alignment count
            gaps_seq = sequence.count("-")
            gaps_alignment += gaps_seq
            
        # Calculate the average number of gaps per sequence
        avg_gaps = round(gaps_alignment / nr_sequences, 2)

        # Add the counts for the alignment to the dataframe using pd.concat
        alignment_name = alignment_file.split("/")[-1] # Get the name of the alignment file
        # Remove .fasta extension if present
        if alignment_name.endswith(".fasta"):
            alignment_name = alignment_name[:-6]
        alignment_results_df = pd.DataFrame({"alignment": [alignment_name], "nr_sequences": nr_sequences, "sequence_length": sequence_length, "sp_score": sp_score, "avg_gaps": avg_gaps, "stop_codons": [stop_codons_total], "incomplete_codons_ends": [incomplete_codons_ends_total], "incomplete_codons_internal": [incomplete_codons_internal_total]})
        results_df = pd.concat([results_df, alignment_results_df], ignore_index=True)
    
    # Write the results to a csv file
    results_df.to_csv(output_file, index=False)
    
    
    ### Make plot with 2x3 bargraphs, one for each metric
    results_df["alignment"] = results_df["alignment"].str.split("_").str[0]
    results_df = results_df.drop(columns=["sequence_length"])
    
    # Add color column 
    results_df["color"] = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7"]
    
    # Rename columns
    results_df = results_df.rename(columns={
        "nr_sequences": "Nr. of Sequences",
        "sp_score": "Sum-of-Pairs Score",
        "avg_gaps": "Average Gaps per Sequence",
        "stop_codons": "Total Stop Codons",
        "incomplete_codons_ends": "Total Incomplete Codons Terminal",
        "incomplete_codons_internal": "Total Incomplete Codons Internal"
    })

    # Set font to Arial and fontsize to 14
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["font.size"] = 12

    # Set up subplots
    fig, ax = plt.subplots(2, 3, figsize=(12, 8))
    
    # Flatten 2D array of axes for easier iteration
    ax = ax.ravel()
    
    # Letters for subplots
    letters = ["A", "B", "C", "D", "E", "F"]
    
    def title_to_sentence_case(text):
        return text[0].upper() + text[1:].lower()
    
    for i, metric in enumerate(["Nr. of Sequences", "Sum-of-Pairs Score", "Average Gaps per Sequence", "Total Stop Codons", "Total Incomplete Codons Terminal", "Total Incomplete Codons Internal"]):
        ax[i].bar(results_df["alignment"], results_df[metric], color=results_df["color"])
        ax[i].set_title(metric, fontweight="bold")
        ax[i].set_ylabel(title_to_sentence_case(metric))
        ax[i].tick_params(axis="x", rotation=45) # Rotate x labels
        
        # Add letters A-F to the top-left corner of each subplot
        ax[i].text(0.03, 0.97, letters[i], transform=ax[i].transAxes, fontsize=16, fontweight="bold", va="top", ha="left")
        
        # Add values on top of each bar (in exponential notation for SP score)
        for j, value in enumerate(results_df[metric]):
            if metric == "Sum-of-Pairs Score":
                # Format the value in exponential notation and split into two lines
                exp_notation = f"{value:.2e}"
                base, exponent = exp_notation.split("e")
                formatted_value = f"{base}\ne{exponent}"  # Separate into two lines
                ax[i].text(j, value, formatted_value, ha="center", va="bottom", fontsize=12)  # Adjust fontsize if needed
            else:
                ax[i].text(j, value, str(value), ha="center", va="bottom")
                
            # Extend the y-axis for values on top of bars and subplot labels
            if metric == "Sum-of-Pairs Score":
                ax[i].set_ylim(0, results_df[metric].max() * 1.4)  # Extend more for this metric
            else:
                ax[i].set_ylim(0, results_df[metric].max() * 1.3)
        
    plt.tight_layout()
    plt.savefig(args.plot)
    
        
