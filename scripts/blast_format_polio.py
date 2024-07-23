import pandas as pd
import re

## Read and process the blast output and sequence metadata

blast_results = pd.read_csv("../data/blast_polio_results.csv", names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov"])
sequence_metadata = pd.read_csv("../data/polio_metadata.csv")

# Select relevant columns in sequence metadata
sequence_metadata = sequence_metadata[["GenBank Accessions", "Genome Name",  "Genome ID", "Genome Status", "Strain", "Completion Date", "Sequencing Platform", "Size", "Isolation Source", "Collection Date", "Isolation Country", "Geographic Group", "Geographic Location", "Host Name", "Lab Host"]]

# Filter sequence metadata to exclude samples from labs / cell culture
sequence_metadata = sequence_metadata[sequence_metadata["Lab Host"].isna()]

# Remove "accn|" string in blast_results qseqid column and rename it
blast_results["qseqid"] = blast_results["qseqid"].str.removeprefix("accn|")
blast_results = blast_results.rename(columns={"qseqid": "GenBank Accessions"})

# Inner merge the two dataframes
blast_results = blast_results.merge(sequence_metadata, on="GenBank Accessions", how="inner")

# Sort by pident
blast_results = blast_results.sort_values(by="pident", ascending=False)

# Filter out query sequences shorter than 500 bp
blast_results = blast_results[blast_results["Size"] >= 500]

# Filter by match length (try 300)
blast_results = blast_results[blast_results["length"] >= 300]

# Extract serotype from genome name
def extract_serotype(genome_name):
    match_polio = re.search("poliovirus [123]", genome_name, re.IGNORECASE)
    if bool(match_polio):
        return match_polio.group().title().replace(" ", "")
    else: 
        return "Unidentified Poliovirus"

blast_results["Serotype"] = blast_results["Genome Name"].map(extract_serotype)

# Count the number of occurences of each serotype
blast_results["Serotype"].value_counts()

## Map sequences identified by BLAST to the corresponding gene (by sequence location)
# Create dataframe of gene locations in the Sabin genomes from the .gb files

# Get list of EV-C .gbff files
import os
evc_files = os.listdir("../data/evc_reference_genomes") 
evc_files = [file for file in evc_files if ".gbff" in file] # Filter for files with .gbff extension

# Get paths of .gbff files
evc_paths = ["../data/evc_reference_genomes/" + file for file in evc_files]

# Create dictionary
evc_files = [file.replace(".gbff", "") for file in evc_files] # Remove .gb extension
evc_dict = dict(zip(evc_files, evc_paths))

# Order dict alphabetically
evc_dict = dict(sorted(evc_dict.items()))


# Extract gene locations from .gbff files using BioPython
# Output: dataframe with EV-C serotype, gene name, gene start, gene end

from Bio import SeqIO

gene_loc_df = pd.DataFrame(columns=["evc_serotype", "gene", "start", "end"])
for evc_serotype in evc_dict.keys():
    gb_path = evc_dict[evc_serotype]
    for record in SeqIO.parse(gb_path, "genbank"):
        for feature in record.features:
            if feature.type == "mat_peptide":
                # Use pd.concat to append rows to the dataframe
                gene_loc_df = pd.concat([gene_loc_df, pd.DataFrame({"evc_serotype": [evc_serotype], "gene": [feature.qualifiers["product"][0]], "start": [feature.location.start], "end": [feature.location.end]})], ignore_index=True)      
            elif feature.type == "5'UTR":
                gene_loc_df = pd.concat([gene_loc_df, pd.DataFrame({"evc_serotype": [evc_serotype], "gene": ["5'UTR"], "start": [feature.location.start], "end": [feature.location.end]})], ignore_index=True)
            elif feature.type == "3'UTR":
                gene_loc_df = pd.concat([gene_loc_df, pd.DataFrame({"evc_serotype": [evc_serotype], "gene": ["3'UTR"], "start": [feature.location.start], "end": [feature.location.end]})], ignore_index=True)
                

# Add column to dataframe: translate sseqid to EV-C reference serotypes
# Create dictionary to map GenBank accession to EV-C serotype
evc_translation_dict = {}
for path in evc_dict.values():
    serotype = path.split("/")[-1].replace(".gbff", "")
    for record in SeqIO.parse(path, "genbank"):
        evc_translation_dict[record.id] = serotype
        
# Map to sseqid
blast_results["Reference EV-C Serotype"] = blast_results["sseqid"].map(evc_translation_dict)
    
    
## Add two columns to blast_results: sgene_start and sgene_end
# Find genes in which the start and end positions of the subject sequence fall
# If both start and end positions lie in the same gene, assign that gene to both sgene_start and sgene_span

def map_gene_start(row):
    evc_reference = row["Reference EV-C Serotype"]
    start = row["sstart"]
    try:
        gene_start = gene_loc_df[(gene_loc_df["evc_serotype"] == evc_reference) & (gene_loc_df["start"] <= start) & (gene_loc_df["end"] >= start)]["gene"].values[0]
    except:
        gene_start = "NA"
    return gene_start

def map_gene_end(row):
    evc_reference = row["Reference EV-C Serotype"]
    end = row["send"]
    try:
        gene_end = gene_loc_df[(gene_loc_df["evc_serotype"] == evc_reference) & (gene_loc_df["start"] <= end) & (gene_loc_df["end"] >= end)]["gene"].values[0]
    # Problem: running into stop codons that are not assigned to a gene
    # Return NA if no gene is found
    except:
        try:
            # Try if gene sequence -3 bp leads to a result -> in that case it's the stop codon
            gene_end = gene_loc_df[(gene_loc_df["evc_serotype"] == evc_reference) & (gene_loc_df["start"] <= end-3) & (gene_loc_df["end"] >= end-3)]["gene"].values[0]
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
    "Reference EV-C Serotype",
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

# Kick out the EV-C116 matches missing gene annotations for the start OR end gene
blast_results = blast_results[blast_results["sgene_start"] != "NA"]
blast_results = blast_results[blast_results["sgene_end"] != "NA"]

# Read dataframe of percentage identity between EV-C reference genomes (genes) and Sabin strains
pident_df = pd.read_csv("../data/pident_evc_sabin.csv")


### Try to calculate sequence identity for the specific sequence fragments instead of on rough gene level?
# Problem: How to find correct gene start and end positions in the EV-C reference genomes?
# (1.) Possible solution: Find best match  in the EV-C reference genome
# (1.) Possible problems: Programatically complicated (and perhaps computationally expensive?), might not find correct match for short and dissimilar sequences
# (2.) Possible solution: Calculate sequence identity (as an average weighted by sequence length) for all the genes that the sequence fragment overlaps with

test_df = blast_results.head(70)

def get_weighted_pident(row):
    evc_genome_order = ["5'UTR", "VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "3'UTR"]
    evc_ref = row["Reference EV-C Serotype"]
    polio_ref = row["Serotype"].replace("Poliovirus", "Sabin")
    subset_df = pident_df[(pident_df["evc_ref"] == evc_ref) & (pident_df["sabin_strain"] == polio_ref)]
    gene_start = row["sgene_start"]
    gene_end = row["sgene_end"]
    genes = evc_genome_order[evc_genome_order.index(gene_start):evc_genome_order.index(gene_end)+1]
    
    # Calculate sum of gene lengths
    gene_length_sum = subset_df[subset_df["gene"].isin(genes)]["evc_gene_length"].sum()
    
    # Calculate weighted average for pident
    weighted_pident = 0
    for gene in genes:
        pident = subset_df.loc[subset_df["gene"] == gene]["pident_ref"].values[0]
        gene_length = subset_df.loc[subset_df["gene"] == gene]["sabin_gene_length"].values[0]
        weighted_pident += pident * gene_length / gene_length_sum
        
    return weighted_pident  

# Remove unidentified polioviruses for now, as these have no reference Sabin strain to calculate expected sequence identities for
blast_results = blast_results[blast_results["Serotype"] != "Unidentified Poliovirus"]

# Run the function on the blast_results dataframe
blast_results["weighted_pident"] = blast_results.apply(get_weighted_pident, axis=1)

# Calculate pident difference between observed and expected
blast_results["weighted_pident"] = blast_results["weighted_pident"]*100 # Convert to percentage
blast_results["pident_diff"] = blast_results["pident"] - blast_results["weighted_pident"]

# Sort by pident_diff
blast_results = blast_results.sort_values(by="pident_diff", ascending=False)

# Filter for pident values at least 5 percentage points bigger than pident_ref and absolute pident at least 85%
blast_filtered_polio = blast_results[(blast_results["pident_diff"] >= 5) & (blast_results["pident"] >= 85)]