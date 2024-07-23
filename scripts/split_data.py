### This script splits the metadata and sequences into poliovirus and non-poliovirus EV-C files.

import pandas as pd
from Bio import SeqIO
import argparse

if __name__ == "__main__":
    
    # Parse arguments
    argparser = argparse.ArgumentParser(description="Split metadata and sequences into poliovirus and non-poliovirus EV-C.")
    argparser.add_argument("--metadata", help="Metadata file in csv format.")
    argparser.add_argument("--sequences", help="Sequences file in fasta format.")
    argparser.add_argument("--genotyping", help="RIVM Genotyping Tool results file in csv format.")
    argparser.add_argument("--output_polio_metadata", help="Output file for the poliovirus metadata.")
    argparser.add_argument("--output_polio_sequences", help="Output file for the poliovirus sequences.")
    argparser.add_argument("--output_npevc_metadata", help="Output file for the non-poliovirus EV-C metadata.")
    argparser.add_argument("--output_npevc_sequences", help="Output file for the non-poliovirus EV-C sequences.")
    
    args = argparser.parse_args()

    ## Read and process RIVM Genotyping Tool results
    genotyping = pd.read_csv(args.genotyping)
    genotyping.rename(columns={"name": "Accession"}, inplace=True) # Rename "name" column to "Accession" to match metadata
    
    # Convert "CV-" to "Coxsackievirus ", "PV-" to "Poliovirus ", and "EV-" to "Enterovirus " in the "VP1 type" column
    genotyping["VP1 type"] = genotyping["VP1 type"].str.replace("CV-", "Coxsackievirus ", regex=False)
    genotyping["VP1 type"] = genotyping["VP1 type"].str.replace("PV-", "Poliovirus ", regex=False)
    genotyping["VP1 type"] = genotyping["VP1 type"].str.replace("EV-", "Enterovirus ", regex=False)
    genotyping["VP1 type"] = genotyping["VP1 type"].str.replace("A24v", "A24", regex=False)
    
    # If inferred VP1 type ist not an EV-C serotype, replace with NA
    evc_serotypes = [
        "Coxsackievirus A1", 
        "Coxsackievirus A11",
        "Coxsackievirus A13",
        "Coxsackievirus A15",
        "Coxsackievirus A17",
        "Coxsackievirus A18",
        "Coxsackievirus A19",
        "Coxsackievirus A20",
        "Coxsackievirus A21",
        "Coxsackievirus A22",
        "Coxsackievirus A24",
        "Enterovirus C95",
        "Enterovirus C96",
        "Enterovirus C99",
        "Enterovirus C102",
        "Enterovirus C104",
        "Enterovirus C105",
        "Enterovirus C109",
        "Enterovirus C113",
        "Enterovirus C117",
        "Enterovirus C118",
        "Poliovirus 1",
        "Poliovirus 2",
        "Poliovirus 3"
    ]
    
    # Set VP1 type to NA if not in EV-C serotypes
    genotyping["VP1 type"] = genotyping["VP1 type"].apply(
        lambda x: x if pd.notna(x) and any([serotype in x for serotype in evc_serotypes]) else pd.NA
    )
        
    # Incorporate the additional VP1 types infered by the RIVM Genotyping Tool into the metadata (use "VP1 types" column)
    metadata = pd.read_csv(args.metadata)
    
    metadata = metadata.merge(genotyping[["Accession", "VP1 type"]], on="Accession", how="left")
    
    
    # If serotype contains "Unidentified" and VP1 type is not NA, replace serotype with VP1 type
    metadata["serotype"] = metadata.apply(
        lambda x: x["VP1 type"] if "Unidentified" in x["serotype"] and not pd.isna(x["VP1 type"]) else x["serotype"],
        axis=1
    )


    ## Split metadata and sequences into poliovirus and non-poliovirus EV-C

    # Get polio metadata by filtering for serotypes that contain "Poliovirus" or "Pv"
    metadata_polio = metadata[metadata["serotype"].str.contains("Poliovirus|Pv")]

    # Get non-polio EV-C metadata by filtering for serotypes that are NOT in polio metadata
    metadata_npevc = metadata[~metadata["serotype"].str.contains("Poliovirus|Pv")]

    # Check if sum of polio and non-polio metadata equals total sequences in metadata
    print(
        "Number of sequences in metadata (after initial filtering): ", len(metadata),
        "\nNumber of polio sequences: ", len(metadata_polio),
        "\nNumber of non-polio EV-C sequences: ", len(metadata_npevc)
    )

    # Read sequences
    sequences = list(SeqIO.parse(args.sequences, "fasta"))

    # Split sequences into two files based on metadata
    sequences_polio = []
    sequences_npevc = []
    for seq in sequences:
        if seq.id in metadata_polio["Accession"].to_list(): 
            sequences_polio.append(seq)
        elif seq.id in metadata_npevc["Accession"].to_list():
            sequences_npevc.append(seq)
            
    # Write sequences to files
    SeqIO.write(sequences_polio, args.output_polio_sequences, "fasta")
    SeqIO.write(sequences_npevc, args.output_npevc_sequences, "fasta")

    # Write metadata to files
    metadata_polio.to_csv(args.output_polio_metadata, index=False)
    metadata_npevc.to_csv(args.output_npevc_metadata, index=False)