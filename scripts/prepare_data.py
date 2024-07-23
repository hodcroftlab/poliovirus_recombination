### Prepare the data downloaded from NCBI Virus for the analysis: read in the metadata and extract the serotype information to split the metadata as well as sequence data (fasta) into poliovirus and non-poliovirus EV-C files
# Input: all_evc_metadata.csv, all_evc_sequences.fasta
# Output: polio_metadata.csv, polio_sequences.fasta, npevc_metadata.csv, npevc_sequences.fasta

import pandas as pd
import country_converter as coco
import re
from Bio import SeqIO
import json
import argparse

if __name__ == "__main__":

    # Parse arguments
    argparser = argparse.ArgumentParser(description="Prepare the data downloaded from NCBI Virus for the analysis: read in the metadata and extract the serotype information to split the metadata as well as sequence data (fasta) into poliovirus and non-poliovirus EV-C files.")
    argparser.add_argument("--metadata", help="Metadata file in csv format.")
    argparser.add_argument("--sequences", help="Sequences file in fasta format.")
    argparser.add_argument("--serotype_corrections", help="A JSON file with manual serotype corrections.")
    argparser.add_argument("--output_metadata", help="Output file for the filtered metadata.")
    argparser.add_argument("--output_sequences", help="Output file for the filtered sequences.")
    argparser.add_argument("--output_unidentified_metadata", help="Output file for the metadata file with unidentified serotypes.")
    argparser.add_argument("--output_unidentified_sequences", help="Output file for the sequence file with unidentified serotypes.")
    argparser.add_argument("--min_length", help="Minimum length of sequences to keep.", type=int, required=False)
    
    args = argparser.parse_args()
    
    ## Read and filter metadata

    raw_metadata = pd.read_csv(args.metadata)
    
    # Filter step 1: filter out sequences from labs / cell culture
    raw_metadata["Isolation_Source"].value_counts() # Check the number of samples from labs / cell culture in the data
    metadata = raw_metadata[raw_metadata["Isolation_Source"] != "cell"]

    # Filter step 2: filter out samples with certain keywords in GenBank title
    keywords = ["unverified", "synthetic", "patent", "cold-adapted", "defective", "identification", "microarray", "oligonucleotide", "pharmaceutical", "modified", "useful", "nonfunctional", "JP 2000503551", "KR 1020140145575", "JP 2015091247", "KR 1020150046348", "JP 2015528285", "JP 2011204261", "JP 2018191645"]
    metadata = metadata[~metadata["GenBank_Title"].str.contains("|".join(keywords), case=False)] # Case-insensitive search

    # Filter step 3: filter out very short sequences (less than min_length)
    # First, filter sequences based on steps 1 and 2
    sequences = list(SeqIO.parse(args.sequences, "fasta")) # Read all sequences
    sequences_filtered = []
    for seq in sequences:
        if seq.id in metadata["Accession"].to_list(): 
            sequences_filtered.append(seq)
            
    # Filter sequences based on length
    min_length = args.min_length if args.min_length else 0
    sequences_filtered = [seq for seq in sequences_filtered if len(seq.seq) >= min_length]
    
    # Write filtered sequences to file
    SeqIO.write(sequences_filtered, args.output_sequences, "fasta")
    
    # Re-adjust metadata accordingly
    metadata = metadata[metadata["Accession"].isin([seq.id for seq in sequences_filtered])]
    
    ## Add continent column to metadata using country_converter
    metadata["continent"] = coco.convert(names=metadata["Country"].tolist(), to="continent_7")
    metadata["continent"] = metadata["continent"].replace("not found", pd.NA) # Convert "not found" to NA

    ## Extract serotypes from metadata

    # Function to extract serotype from organism name, genotype and/or title
    def extract_serotype(genome_name, genotype, gb_title):
        
        # Ensure genome_name and genotype are strings
        genome_name = str(genome_name) if genome_name else ""
        genotype = str(genotype) if genotype else ""
        
        # Check for Coxsackievirus A19/22 first as it's more specific
        match = re.search(r"Coxsackievirus A19/22", genome_name, re.IGNORECASE)
        if match:
            return match.group().title()
        
        # Check for Coxsackievirus A followed by a number
        match = re.search(r"Coxsackievirus A\d+", genome_name, re.IGNORECASE)
        if match:
            return match.group().title()
        
        # Check for Enterovirus C followed by a number
        match = re.search(r"Enterovirus C\d+", genome_name, re.IGNORECASE)
        if match:
            return match.group().title()
        
        # Check for Poliovirus followed by a specific number
        match = re.search(r"poliovirus [123]", genome_name, re.IGNORECASE)
        if match:
            return match.group().title()
        
        # Check for wild poliovirus type 3 specifically
        if re.search(r"wild poliovirus type 3", genome_name, re.IGNORECASE):
            return "Poliovirus 3"
        
        # Check for general poliovirus
        if re.search(r"poliovirus", genome_name, re.IGNORECASE):
            # Check for genotype patterns
            if re.search(r"recombinant", genotype, re.IGNORECASE):
                return genotype.title()
            return "Unidentified Poliovirus"
        
        # Check for general enterovirus and additional genotype info
        if re.search(r"enterovirus", genome_name, re.IGNORECASE):
            # Check for genotype patterns
            match = re.search(r"A\d+", genotype, re.IGNORECASE)
            if match:
                return "Coxsackievirus " + match.group().title()
            
            match = re.search(r"C\d+", genotype, re.IGNORECASE)
            if match:
                return "Enterovirus " + match.group().title()
            
            match = re.search(r"CAV\d+", genotype, re.IGNORECASE)
            if match:
                return match.group().replace("CAV", "Coxsackievirus A")
            
            match = re.search(r"PV\d+", genotype, re.IGNORECASE)
            if match:
                return match.group().replace("PV", "Poliovirus ")
            
            match = re.search(r"poliovirus \d+", genotype, re.IGNORECASE)
            if match:
                return match.group().title()
            
            # Also check GenBank title column for poliovirus and enterovirus C
            match = re.search(r"PV-\d+", gb_title, re.IGNORECASE)
            if match:
                return match.group().replace("PV-", "Poliovirus ")
            
            match = re.search(r"PV\d+", gb_title, re.IGNORECASE)
            if match:
                return match.group().replace("PV", "Poliovirus ")
            
            if re.search(r"PV|polio", gb_title, re.IGNORECASE):
                return "Unidentified Poliovirus"
            
            match = re.search(r"EVC/C\d+", gb_title, re.IGNORECASE)
            if match:
                return match.group().replace("EVC/", "Enterovirus ")
            
            match = re.search(r"CSV/A\d+", gb_title, re.IGNORECASE)
            if match:
                return match.group().replace("CSV/", "Coxsackievirus ")
            
            match = re.search(r"CVA-\d+", gb_title, re.IGNORECASE)
            if match:
                return match.group().replace("CVA-", "Coxsackievirus A")
            
            return "Unidentified Enterovirus C"
        
        return "NA"

    # Apply the function to the dataframe
    metadata["serotype"] = metadata.apply(lambda x: extract_serotype(x["Organism_Name"], x["Genotype"], x["GenBank_Title"]), axis=1)

    # Load manually created serotype corrections and map to dataframe
    serotype_corrections = json.loads(open(args.serotype_corrections).read())
    metadata["serotype"] = metadata["Accession"].map(serotype_corrections).fillna(metadata["serotype"])

    # ON604XXX: VDPV, but unclear which type --> Unidentified Poliovirus
    metadata.loc[metadata["Accession"].str.contains("ON604"), "serotype"] = "Unidentified Poliovirus"

    # Count occurence of each serotype
    print(metadata["serotype"].value_counts())
    
    
    ## Save filtered metadata to file
    metadata.to_csv(args.output_metadata, index=False)
    
    
    ## Save unidentified polio and non-polio EV-C sequences to a separate file
    metadata_unidentified = metadata[metadata["serotype"].str.contains("Unidentified")]
    sequences_unidentified = []
    for seq in sequences_filtered:
        if seq.id in metadata_unidentified["Accession"].to_list(): 
            sequences_unidentified.append(seq)

    
    SeqIO.write(sequences_unidentified, args.output_unidentified_sequences, "fasta")
    metadata_unidentified.to_csv(args.output_unidentified_metadata, index=False)
    
    
    ## Print message giving the number of sequences before and after initial filtering
    print(
        "\nNumber of sequences in original data: ", len(raw_metadata),
        "\nNumber of sequences after filtering: ", len(metadata),
        "\nOf which have an unidentified serotypes: ", len(metadata_unidentified)
    )