### This script assigns recombination status of each sequence to the metadata
### Recombination status is assigned based on sorting of SimPlot results

import pandas as pd
import os
import argparse

argparser = argparse.ArgumentParser(description="Assign recombination status of each sequence to metadata based on sorted SimPlots.")
argparser.add_argument("--metadata_in", help="Input CSV file containing metadata.")
argparser.add_argument("--metadata_out", help="Output CSV file containing updated metadata with recombination status.")
args = argparser.parse_args()

# Read metadata
metadata = pd.read_csv(args.metadata_in)

# Get sequence IDs corresponding to recombination status from folder names
recombination_status_df = pd.DataFrame(columns=["accession", "recombination_status"])
for folder in os.listdir("plots/simplots/all_sorted"):
    accessions = os.listdir(f"plots/simplots/all_sorted/{folder}")
    accessions = [accn.split("_")[1] for accn in accessions] # Extract the sequence IDs
    
    # Create a DataFrame for this folder
    folder_df = pd.DataFrame({
        "accession": accessions, 
        "recombination_status": [folder] * len(accessions)
    })
    # Concatenate to the main DataFrame
    recombination_status_df = pd.concat([recombination_status_df, folder_df], ignore_index=True)

# Merge metadata with recombination status
metadata = metadata.rename(columns={"Accession": "accession"})
updated_metadata = pd.merge(metadata, recombination_status_df, on="accession", how="left")

# Print sequence IDs with missing recombination status
missing_recombination_status = updated_metadata[updated_metadata["recombination_status"].isna()]
if not missing_recombination_status.empty:
    print("The following sequences are missing recombination status:")
    print(missing_recombination_status["accession"])
    
    # Assign "non_recombinant" to these sequences (PV1 Mahoney reference sequence + Sabin vaccine strains)
    updated_metadata.loc[updated_metadata["recombination_status"].isna(), "recombination_status"] = "non_recombinant"
    print("Assigned 'non_recombinant' to these sequences.")

# Save as CSV
updated_metadata.to_csv(args.metadata_out, index=False)