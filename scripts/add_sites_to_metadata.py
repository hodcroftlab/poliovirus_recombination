### This script adds the nucleotide status at genetic sites of particular interest (e.g., attenuation determinants, CRE) to the metadata file.

import pandas as pd
import argparse
from Bio import SeqIO

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Add the nucleotide status at genetic sites of particular interest to the metadata file.")
    parser.add_argument("--metadata", help="Metadata file in .csv format.")
    parser.add_argument("--alignment", help="Multiple sequence alignment in fasta format.")
    parser.add_argument("--sites", help="File containing the genetic sites of interest in .csv format (expected columns: nuc_site, feature).")
    parser.add_argument("--output", help="Output directory for the updated metadata file.")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    sites = pd.read_csv(args.sites, sep="\t")
    alignment = list(SeqIO.parse(args.alignment, "fasta"))
    
    metadata = pd.read_csv("../data/metadata/evc_full_genomes_metadata_rivm.csv")
    sites = pd.read_csv("../data/config/sites_of_interest.tsv", sep="\t")
    alignment = list(SeqIO.parse("../data/sequences/alignments/full_genome/assembled/fullgenome_assembled.fasta", "fasta"))
    
    # Convert sites to dictionary: {feature: [nuc_sites]}
    sites_dict = {feature: group["nuc_site"].tolist() for feature, group in sites.groupby("feature")}

    # Initialize list to collect sequence data
    nuc_status_list = []

    # Iterate over each sequence in the alignment
    for record in alignment:
        seq_id = record.id
        seq = str(record.seq)

        entry_dict = {"Accession": seq_id}  # Initialize entry with strain ID
        
        
        ### FIX THE FOLLOWING PART

        # Iterate over the sites of interest
        for feature, nucs in sites_dict.items():
            site_list = []

            for nuc_site in nucs:
                # Ensure the site is within the sequence range
                if nuc_site <= len(seq):
                    nuc = seq[nuc_site - 1]  # Get nucleotide
                    entry = f"{nuc_site}{nuc}"
                else:
                    entry = "NA"

                site_list.append(entry)

        # Add the entry_dict to the list for later conversion
        nuc_status_list.append(entry_dict)

    # Convert the collected list into a DataFrame
    nuc_status_df = pd.DataFrame(nuc_status_list)

    # Merge with metadata
    updated_metadata = pd.merge(metadata, nuc_status_df, on="Accession", how="left")

    # Save updated metadata to file
    updated_metadata.to_csv(args.output, index=False)

    print(f"Updated metadata saved to {args.output}.")