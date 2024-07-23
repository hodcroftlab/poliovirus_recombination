import pandas as pd
import re

## Read and process the blast output and sequence metadata

# Read csv files, add column names to the blast output
blast_results = pd.read_csv("../data/blast_evc_results.csv", names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov"]) 
sequence_metadata = pd.read_csv("../data/non_polio_evc_metadata.csv")

# Select relevant columns in sequence metadata
sequence_metadata = sequence_metadata[["GenBank Accessions", "Genome Name",  "Genome ID", "Genome Status", "Strain", "Completion Date", "Sequencing Platform", "Size", "Isolation Source", "Collection Date", "Isolation Country", "Geographic Group", "Geographic Location", "Host Name", "Lab Host"]]

# Filter sequence metadata to exclude samples from labs / cell culture
sequence_metadata = sequence_metadata[sequence_metadata["Lab Host"].isna()]

# Extract serotype from genome name	    
    
def extract_serotype(genome_name):
    match_coxa = re.search("Coxsackievirus A\d+", genome_name, re.IGNORECASE)
    match_evc = re.search("Enterovirus C\d+", genome_name, re.IGNORECASE)
    if bool(match_coxa):
        return match_coxa.group().title().replace(" ", "")
    elif bool(match_evc):
        return match_evc.group().title().replace(" ", "")
    
sequence_metadata["Serotype"] = sequence_metadata["Genome Name"].map(extract_serotype)

# Filter out viruses without serotype information (unidentified polioviruses)
sequence_metadata = sequence_metadata[sequence_metadata["Serotype"].notna()]

# According to the ICTV website, Coxsackievirus A15 and Coxsackievirus A11 are the same serotype => label A15 as A11
sequence_metadata["Serotype"] = sequence_metadata["Serotype"].replace("CoxsackievirusA15", "CoxsackievirusA11")

# No full genome sequences available for Enterovirus C95 => remove from the dataset
sequence_metadata = sequence_metadata[sequence_metadata["Serotype"] != "EnterovirusC95"]

# Count the number of occurences of each serotype
sequence_metadata["Serotype"].value_counts()

# Remove "accn|" string in blast_results qseqid column and rename it
blast_results["qseqid"] = blast_results["qseqid"].str.removeprefix("accn|")
blast_results = blast_results.rename(columns={"qseqid": "GenBank Accessions"})

# "Translate" sseqid to Sabin 1/2/3
sabin_translation_dict = {"AY184219.1": "Sabin1", "AY184220.1": "Sabin2", "AY184221.1": "Sabin3"}

blast_results["Reference Sabin Strain"] = blast_results["sseqid"].map(sabin_translation_dict)
    
# Inner merge the two dataframes
blast_results = blast_results.merge(sequence_metadata, on="GenBank Accessions", how="inner")

# Sort by pident
blast_results = blast_results.sort_values(by="pident", ascending=False)

# Filter out query sequences shorter than 500 bp
blast_results = blast_results[blast_results["Size"] >= 500]

# Filter by match length (try 300)
blast_results = blast_results[blast_results["length"] >= 300]


## Map sequences identified by BLAST to the corresponding gene (by sequence location)
# Create dataframe of gene locations in the Sabin genomes from the .gb files

# Get list of Sabin .gb files
import os
sabin_files = os.listdir("../data/sabin123_reference_genomes") 
sabin_files = [file for file in sabin_files if ".gb" in file] # Filter for files with .gb extension

# Get paths of .gb files
sabin_paths = ["../data/sabin123_reference_genomes/" + file for file in sabin_files]

# Create dictionary
sabin_files = [file.replace(".gb", "") for file in sabin_files] # Remove .gb extension
sabin_dict = dict(zip(sabin_files, sabin_paths))

# Order dict to start with Sabin1
sabin_dict = dict(sorted(sabin_dict.items()))

# Extract gene locations from .gb files using BioPython
# Output: dataframe with Sabin strain, gene name, gene start, gene end

from Bio import SeqIO

gene_loc_df = pd.DataFrame(columns=["sabin_strain", "gene", "start", "end"])
for sabin_strain in sabin_dict.keys():
    gb_path = sabin_dict[sabin_strain]
    for record in SeqIO.parse(gb_path, "genbank"):
        for feature in record.features:
            if feature.type == "mat_peptide":
                # Use pd.concat to append rows to the dataframe
                gene_loc_df = pd.concat([gene_loc_df, pd.DataFrame({"sabin_strain": [sabin_strain], "gene": [feature.qualifiers["product"][0]], "start": [feature.location.start], "end": [feature.location.end]})], ignore_index=True)      
            elif feature.type == "5'UTR":
                gene_loc_df = pd.concat([gene_loc_df, pd.DataFrame({"sabin_strain": [sabin_strain], "gene": ["5'UTR"], "start": [feature.location.start], "end": [feature.location.end]})], ignore_index=True)
            elif feature.type == "3'UTR":
                gene_loc_df = pd.concat([gene_loc_df, pd.DataFrame({"sabin_strain": [sabin_strain], "gene": ["3'UTR"], "start": [feature.location.start], "end": [feature.location.end]})], ignore_index=True)
    
## Add two columns to blast_results: sgene_start and sgene_end
# Find gene in which the start and end positions of the subject sequence fall
# If both start and end positions lie in the same gene, assign that gene to sgene_start and sgene_span
    
def map_gene_start(row):
    serotype = row["Reference Sabin Strain"]
    start = row["sstart"]
    try:
        gene_start = gene_loc_df[(gene_loc_df["sabin_strain"] == serotype) & (gene_loc_df["start"] <= start) & (gene_loc_df["end"] >= start)]["gene"].values[0]
    except:
        gene_start = "NA"
    return gene_start

def map_gene_end(row):
    serotype = row["Reference Sabin Strain"]
    end = row["send"]
    try:
        gene_end = gene_loc_df[(gene_loc_df["sabin_strain"] == serotype) & (gene_loc_df["start"] <= end) & (gene_loc_df["end"] >= end)]["gene"].values[0]
    except: 
        gene_end = "NA"
    return gene_end
    
# Run the function on the blast_results dataframe
blast_results["sgene_start"] = blast_results.apply(map_gene_start, axis=1)
blast_results["sgene_end"] = blast_results.apply(map_gene_end, axis=1)

# Reorder columns in blast_results
blast_results = blast_results[[
    "GenBank Accessions",
    "Genome Name",
    "Serotype",
    "Genome Status",
    "Reference Sabin Strain",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sgene_start",
    "sgene_end",
    "evalue",
    "bitscore",
    "qcov",
    "Host Name",
    "Size",
    "Isolation Source",
    "Collection Date",
    "Isolation Country",
    "Geographic Group",
    "Geographic Location",
    "Genome ID",
    "Strain",
    "Completion Date",
    "Sequencing Platform",
    "sseqid"]]



# Read dataframe of percentage identity between EV-C reference genomes (genes) and Sabin strains
pident_df = pd.read_csv("../data/pident_evc_sabin.csv")


### Try to calculate sequence identity for the specific sequence fragments instead of on rough gene level?
# Problem: How to find correct gene start and end positions in the EV-C reference genomes?
# (1.) Possible solution: Find best match  in the EV-C reference genome
# (1.) Possible problems: Programatically complicated (and perhaps computationally expensive?), might not find correct match for short and dissimilar sequences
# (2.) Possible solution: Calculate sequence identity (as an average weighted by sequence length) for all the genes that the sequence fragment overlaps with


def get_weighted_pident(row):    
    evc_genome_order = ["5'UTR", "VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "3'UTR"]
    evc_ref = row["Serotype"]
    sabin_strain = row["Reference Sabin Strain"]
    gene_start = row["sgene_start"]
    gene_end = row["sgene_end"]
    genes = evc_genome_order[evc_genome_order.index(gene_start):evc_genome_order.index(gene_end)+1]
    
    # Calculate sum of gene lengths
    gene_length_sum = 0
    for gene in genes:
        gene_length_sum += pident_df[(pident_df["evc_ref"] == evc_ref) & (pident_df["gene"] == gene) & (pident_df["sabin_strain"] == sabin_strain)]["evc_gene_length"].values[0]
    
    # Calculate weighted average for pident
    weighted_pident = 0
    for gene in genes:
        pident = pident_df[(pident_df["evc_ref"] == evc_ref) & (pident_df["gene"] == gene) & (pident_df["sabin_strain"] == sabin_strain)]["pident_ref"].values[0]
        gene_length = pident_df[(pident_df["evc_ref"] == evc_ref) & (pident_df["gene"] == gene) & (pident_df["sabin_strain"] == sabin_strain)]["evc_gene_length"].values[0]
        weighted_pident += pident * gene_length / gene_length_sum
        
    return weighted_pident

# Run the function on the blast_results dataframe
blast_results["weighted_pident"] = blast_results.apply(get_weighted_pident, axis=1)

# Convert weighted_pident column to percentage
blast_results["weighted_pident"] = blast_results["weighted_pident"]*100

# Add new column: pident_diff
blast_results["pident_diff"] = blast_results["pident"] - blast_results["weighted_pident"]

# Filter for pident values at least 5 percentage points bigger than pident_ref and absolute pident at least 85%
blast_filtered_npevc = blast_results[(blast_results["pident_diff"] >= 5) & (blast_results["pident"] >= 85)]

