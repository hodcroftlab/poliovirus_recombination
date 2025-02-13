import pandas as pd
import os
import argparse

argparser = argparse.ArgumentParser(description="Identify breakpoints in the genome based on sequence similarity comparison.")
argparser.add_argument("--similarity_data", nargs="+", help="List of CSV files containing pairwise sequence similarity results.")
argparser.add_argument("--highest_similarity_threshold", type=float, default=0.98, help="Sequence similarity threshold above which the highest match is always assigned as parent serotype.")
argparser.add_argument("--high_similarity_threshold", type=float, default=0.95, help="Sequence similarity threshold above which the highest match is assigned as parent serotype if the difference to the next highest similarity is above the low_difference_threshold.")
argparser.add_argument("--low_similarity_threshold", type=float, default=0.88, help="Sequence similarity threshold above which the highest match is assigned as parent serotype if the difference to the next highest similarity is above the high_difference_threshold.")
argparser.add_argument("--high_difference_threshold", type=float, default=0.08, help="Threshold for difference to next best match; used when similarity is > high_similarity_threshold.")
argparser.add_argument("--low_difference_threshold", type=float, default=0.02, help="Threshold for difference to next best match; used when similarity is < high_similarity_threshold and > low_similarity_threshold.")
argparser.add_argument("--output", help="Output CSV file containing breakpoints matched to sequences.")
args = argparser.parse_args()

# Initialize dataframe for breakpoints
breakpoints_df = pd.DataFrame(columns=["sequence", "breakpoint"])

similarity_data = args.similarity_data

highest_similarity_threshold = args.highest_similarity_threshold
high_similarity_threshold = args.high_similarity_threshold
low_similarity_threshold = args.low_similarity_threshold
high_difference_threshold = args.high_difference_threshold
low_difference_threshold = args.low_difference_threshold


# similarity_data = ["../data/recombination_analysis/custom_simplots/results/AY948201.1_similarity_results.csv"]
# high_similarity_threshold = 0.97
# low_similarity_threshold = 0.88
# high_difference_threshold = 0.08
# low_difference_threshold = 0.02

for csv in similarity_data:
    print("Inferring breakpoints for", csv.split("/")[-1].split("_")[0], "...")
    # Read data from the similarity results files
    data = pd.read_csv(csv)

    # Remove non-consensus sequences (identifiable by presence of a dot "." in the seq1 column)
    data = data[~data["seq1"].str.contains("\.")]

    # For each step, find the sequence with the lowest p-distance and calculate the difference to the next lowest p-distance
    # Make new dataframe for this

    data = data.sort_values(["step", "pdist"]).reset_index(drop=True)
    idx = data.groupby("step").idxmax(numeric_only=True)["similarity"]
    winners = data.loc[idx].reset_index(drop=True)

    # Add a column with the difference to the next highest similarity
    winners["diff_to_next"] = winners["similarity"].values - data.loc[idx + 1, "similarity"].values # works because data is sorted by pdist
    
    
    
    ## FIRST STEP OF ASSIGNING PARENTS

    # Make an extra column indicating confidence of the winner: 
    # A) True if similarity > highest_similarity_threshold, irrespective of diff_to_next
    # B) True if similarity > high_similarity_threshold & diff_to_next > low_difference_threshold
    # C) True if similarity > low_similarity_threshold & diff_to_next > high_difference_threshold
    # Assign parent based on this; if confidence is True, assign parent as serotype1, else assign "?" as parent
    winners["confidence"] = (
        (winners["similarity"] >= highest_similarity_threshold) # Very high similarity; always assign as parent
        | ((winners["similarity"] >= high_similarity_threshold) & (winners["diff_to_next"] >= low_difference_threshold)) # High similarity; lower difference is tolerated
        | ((winners["similarity"] >= low_similarity_threshold) & (winners["diff_to_next"] >= high_difference_threshold)) # Low similarity; higher difference needed
    )
    winners["parent"] = ["?" if not conf else serotype for conf, serotype in zip(winners["confidence"], winners["serotype1"])]
    
    # Copy the column to a new column to keep track of the original parent assignment
    winners["original_parent"] = winners["parent"]
    

    
    ## SECOND STEP OF ASSIGNING PARENTS
    
    # Take sequence context into account
    # Find small groups of consecutive serotypes
    # Apply stricter criteria to these; otherwise, classify as noise
    
    min_consecutive = 5  # Threshold to identify small groups

    winners["group"] = (winners["parent"] != winners["parent"].shift()).cumsum()  # Identify groups of consecutive types
    grouped = winners.groupby("group")  # Group by consecutive runs
    
    # Create a column that counts the size of each group
    winners["group_size"] = grouped["parent"].transform("size")
    
    # Identify groups smaller than the threshold
    all_groups = winners["group"].unique()
    small_groups = winners.loc[winners["group_size"] < min_consecutive, "group"].unique()
    
    # First step: sort small groups into "real" and "noise" based on sequence similarity
    noise_groups = []
    real_groups = []
    
    for group in small_groups:
        # Get average similarity of the group
        # Apply stricter criteria to find "real" small groups
        group_sim = winners.loc[winners["group"] == group, "similarity"].mean()
        group_diff = winners.loc[winners["group"] == group, "diff_to_next"].mean()
        if ((group_sim >= high_similarity_threshold) & (group_diff >= high_difference_threshold)) or ((group_sim >= highest_similarity_threshold) & (group_diff >= low_difference_threshold)):
            real_groups.append(group)
        else:
            noise_groups.append(group)
            
    # Second step: iterate over noise_groups and assign parents based on surrounding serotypes
    for group in noise_groups:
        # Check if previous group exists or if it's the first group
        if group - 1 in all_groups:
            prev_group = group - 1
            if prev_group in noise_groups:
                prev_serotype = "?"
            else:
                prev_serotype = winners.loc[winners["group"] == prev_group, "parent"].values[0]
        else: # Start of the genome
            prev_group = None
            
        # Check if next group exists or if it's the last group
        if group + 1 in all_groups:
            next_group = group + 1
            if next_group in noise_groups:
                next_serotype = "?"
            else:
                next_serotype = winners.loc[winners["group"] == next_group, "parent"].values[0]
        else: # End of the genome
            next_group = None
        
        # If previous group is None, assign serotype of next group, and vice versa
        if prev_group == None:
            winners.loc[winners["group"] == group, "parent"] = next_serotype
        elif next_group == None:
            winners.loc[winners["group"] == group, "parent"] = prev_serotype
        else:
            # If both surrounding groups have the same serotype, assign that serotype to the group
            if prev_serotype == next_serotype:
                winners.loc[winners["group"] == group, "parent"] = prev_serotype
            else: # If the surrounding serotypes are different, split the group
                if prev_serotype == "?" or next_serotype == "?": # If either serotype is "?", assign "?" to the group
                    winners.loc[winners["group"] == group, "parent"] = "?"
                else: # If both serotypes are known, split the group in the middle
                    group_size = winners.loc[winners["group"] == group, "group_size"].values[0]
                    # If the group size is even, split in half
                    if group_size % 2 == 0:
                        split_idx = group_size // 2
                        winners.loc[winners["group"] == group, "parent"] = [prev_serotype] * split_idx + [next_serotype] * split_idx
                    else:
                        # If the group size is odd, assign previous serotype to first half and next serotype to second half (excluding the middle step)
                        # Assign the middle step to the serotype with the higher similarity at that step

                        # Find the middle step
                        steps = winners.loc[winners["group"] == group, "step"].values
                        middle_step = steps[len(steps) // 2]
                        
                        # Assign serotypes to the first and second half
                        winners.loc[(winners["group"] == group) & (winners["step"] < middle_step), "parent"] = prev_serotype
                        winners.loc[(winners["group"] == group) & (winners["step"] > middle_step), "parent"] = next_serotype
                        
                        # Assign serotype to the middle step: find the serotype with the highest similarity at the middle step
                        prev_sim = data.loc[(data["step"] == middle_step) & (data["serotype1"] == prev_serotype), "similarity"].values[0]
                        next_sim = data.loc[(data["step"] == middle_step) & (data["serotype1"] == next_serotype), "similarity"].values[0]
                        if prev_sim > next_sim:
                            winners.loc[(winners["group"] == group) & (winners["step"] == middle_step), "parent"] = prev_serotype
                        else:
                            winners.loc[(winners["group"] == group) & (winners["step"] == middle_step), "parent"] = next_serotype
                
                
    ## THIRD STEP OF ASSIGNING PARENTS
    
    # For 1-3 step groups at the very start of the genome, re-assign parents based on next group
    
    # Reassign groups and group sizes based on the updated parent assignments
    winners["group"] = (winners["parent"] != winners["parent"].shift()).cumsum() 
    grouped = winners.groupby("group")  # Group by consecutive runs
    winners["group_size"] = grouped["parent"].transform("size")
    
    # Check if the first group is 1-3 steps long
    first_group = winners["group"].min()
    first_group_size = winners.loc[winners["group"] == first_group, "group_size"].values[0]
    if first_group_size in [1, 2, 3]:
        # Find the next group
        next_group = first_group + 1
        next_serotype = winners.loc[winners["group"] == next_group, "parent"].values[0]
        winners.loc[winners["group"] == first_group, "parent"] = next_serotype

    # Write the results to a CSV file
    path_out = csv.replace("_similarity_results.csv", "_winners.csv").replace("results", "winners")
    winners.to_csv(path_out, index=False)


    ### Assign the breakpoints based on parent column

    # When the parent column changes, assign a breakpoint at the mean position of the two steps

    # Initialize a list to store the breakpoints
    breakpoints = []
    
    # Iterate over the rows in the dataframe
    for i in range(1, len(winners)):
        if winners["parent"][i] != winners["parent"][i-1]:
            # Assign a breakpoint at the mean position of the two steps
            breakpoint = (winners["step"][i] + winners["step"][i-1]) / 2
            breakpoints.append(breakpoint)
    
    # Add to breakpoints dataframe
    sequence = csv.split("/")[-1].split("_")[0]
    breakpoints_df = pd.concat([breakpoints_df, pd.DataFrame({"sequence": [sequence] * len(breakpoints), "breakpoint": breakpoints})])
    
# Write breakpoints to CSV
breakpoints_df = breakpoints_df.sort_values("sequence").reset_index(drop=True)
breakpoints_df.to_csv(args.output, index=False)