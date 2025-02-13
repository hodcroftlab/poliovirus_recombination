### This script corrects unidentified and misidentified serotypes in the metadata file based on the RIVM genotyping results

import pandas as pd
import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Correct unidentified and misidentified serotypes in the metadata file based on the RIVM genotyping results.")
    argparser.add_argument("--metadata", help="Path to CSV/TSV metadata file.", required=True)
    argparser.add_argument("--genotypes", help="Path to CSV RIVM genotyping results file.", required=True)
    argparser.add_argument("--output", help="Output path for corrected metadata file.", required=True)
    args = argparser.parse_args()

    # Read inputs
    rivm_results = pd.read_csv(args.genotypes)
    #rivm_results = pd.read_csv("../data/metadata/rivm_genotyping_results_fullgenome.csv")
    rivm_results = rivm_results[["name", "VP1 type", "VP1 type support"]]
    
    if args.metadata.endswith(".csv"):
        metadata = pd.read_csv(args.metadata)
    elif args.metadata.endswith(".tsv"):
        metadata = pd.read_csv(args.metadata, sep="\t")
    else:
        raise ValueError("Metadata file must be a CSV or TSV file.")
    
    #metadata = pd.read_csv("../data/metadata/evc_full_genomes_metadata.csv")

    # Adjust type names to match metadata
    rivm_results["VP1 type"] = rivm_results["VP1 type"].str.replace("CV-A", "Coxsackievirus A")
    rivm_results["VP1 type"] = rivm_results["VP1 type"].str.replace("EV-C", "Enterovirus C")
    rivm_results["VP1 type"] = rivm_results["VP1 type"].str.replace("PV-", "Poliovirus ")

    # Replace A24v with A24
    rivm_results["VP1 type"] = rivm_results["VP1 type"].str.replace("A24v", "A24")

    # Merge metadata with RIVM results
    rivm_results = rivm_results.rename(columns={"name": "Accession", "VP1 type": "rivm_type", "VP1 type support": "type_support"})
    metadata = metadata.merge(rivm_results, on="Accession", how="left")

    # Filter for rows where serotype is not equal to rivm_type
    mismatches = metadata[metadata["serotype"] != metadata["rivm_type"]]

    # If rivm_type == "Could not assign", keep serotype
    # Else, replace serotype with rivm_type
    mismatches = mismatches[mismatches["rivm_type"] != "Could not assign"]
    mismatches["serotype"] = mismatches["rivm_type"]

    # Replace the changed rows in the metadata
    updated_metadata = pd.concat([mismatches, metadata])   # Concat: changed rows will be at the top
    updated_metadata = updated_metadata.drop_duplicates(subset="Accession")    # Drop duplicates
    
    # Make an extra column with shortened serotype names
    updated_metadata["serotype_short"] = updated_metadata["serotype"].str.replace("Coxsackievirus A", "CVA")
    updated_metadata["serotype_short"] = updated_metadata["serotype_short"].str.replace("Enterovirus C", "EVC")
    updated_metadata["serotype_short"] = updated_metadata["serotype_short"].str.replace("Poliovirus ", "PV")
    
    # Drop rivm_type column
    updated_metadata = updated_metadata.drop("rivm_type", axis=1)

    # Write to file
    if args.output.endswith(".csv"):
        updated_metadata.to_csv(args.output, index=False)
    elif args.output.endswith(".tsv"):
        updated_metadata.to_csv(args.output, sep="\t", index=False)
    else:
        raise ValueError("Output file must be a CSV or TSV file.")