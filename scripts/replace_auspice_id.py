### This script replaces the sequence identifiers in an auspice json file (e.g., from accession ID to GenBank title).

import argparse
import json
import pandas as pd

if __name__ == "__main__":

    # Define argparse arguments
    argparser = argparse.ArgumentParser("Adapt sequence identifiers in an auspice json file: from 'Accession' to 'GenBank Title (Accession)'.")
    argparser.add_argument("--auspice_json", help="Input auspice json file.")
    argparser.add_argument("--metadata", help="Input metadata file.")
    argparser.add_argument("--accn_col", help="Column in metadata file containing accession IDs.")
    argparser.add_argument("--gbtitle_col", help="Column in metadata file containing GenBank titles.")
    argparser.add_argument("--output_auspice_json", help="Output auspice json file.")

    args = argparser.parse_args()

    # Read metadata file
    metadata = pd.read_csv(args.metadata)
    # Add a column with the GenBank title (Accession)
    metadata["new_id"] = metadata[args.gbtitle_col] + " (" + metadata[args.accn_col] + ")"
    # Create a dictionary with the mapping between old and new identifiers
    name_mapping = dict(zip(metadata[args.accn_col], metadata["new_id"]))
        
    # Read auspice json file
    with open(args.auspice_json) as handle:
        auspice_json = json.load(handle)
            
    # Replace sequence identifiers in auspice json file with unique_id using a recursive function
    def replace_name(node):
        if node["name"] in name_mapping:
            node["name"] = name_mapping[node["name"]]
        # Check if the node has children
        if "children" in node:
            for child in node["children"]:
                replace_name(child)

    # Replace the names starting from the root
    replace_name(auspice_json["tree"])
        
    # Write modified auspice json file
    with open(args.output_auspice_json, "w") as handle:
        json.dump(auspice_json, handle)
        

    

