### Part 1: extract nucleotide sequences of genes from genbank files
# Input: gbff files with gene annotations and nucleotide sequences
# Output: Dictionary with gene annotations as keys and nucleotide sequences as values?

## Notes
# Genes/protein products including sequence ranges are annotated as "mat_peptide" in the .gbff files
# 5'-UTR and 3'-UTR are also annotated in the gbff files and need to be included

from Bio import SeqIO
import os

# Create dictionary of .gbff files (remove fasta directory from path)
gb_files = os.listdir("../data/evc_reference_genomes")
gb_files.remove("fasta")
gb_paths = ["../data/evc_reference_genomes/" + file for file in gb_files]

# Remove .gbff extension
gb_files = [file.replace(".gbff", "") for file in gb_files]

# Zip together file names and paths to create dictionary
gb_dict = dict(zip(gb_files, gb_paths))

# Extract gene sequences from .gbff files
def extract_gene_sequences(gb_file):
    gene_dict = {}
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "mat_peptide": # Get protein coding genes
                gene_name = feature.qualifiers["product"][0]
                gene_seq = feature.location.extract(record.seq)
                gene_dict[gene_name] = gene_seq
            elif feature.type == "5'UTR": # Get 5'UTR
                utr5_seq = feature.location.extract(record.seq)
                gene_dict["5'UTR"] = utr5_seq
            elif feature.type == "3'UTR": # Get 3'UTR
                utr3_seq = feature.location.extract(record.seq)
                gene_dict["3'UTR"] = utr3_seq 
            elif feature.type == "source": # Get full genome sequence
                full_seq = feature.location.extract(record.seq)
                if len(full_seq) > 6000: # Check if full sequence is longer than 6000 bp to exclude partial genomes
                    gene_dict["full_genome"] = full_seq
    return gene_dict

# Test the function
extract_gene_sequences(gb_dict["EnterovirusC102"])

# Create dictionary with gene sequences for each .gbff file
evc_gene_dict = {}
for gb_path in gb_dict.values():
    gene_dict = extract_gene_sequences(gb_path)
    evc = [key for key in gb_dict.keys() if gb_dict[key] == gb_path][0]
    evc_gene_dict[evc] = gene_dict
    
# Visually inspect the newly created dictionary
evc_gene_dict
evc_gene_dict.keys()



### Part 2: compute pairwise sequence similarity between EV-C reference genes and Sabin vaccine strains

## Create dictionary of Sabin vaccine strains
# Note: I manually edited the gb files for the Sabin strains to include annotations for the 5' and 3' UTRs

# Get list of Sabin .gb files
sabin_files = os.listdir("../data/sabin123_reference_genomes") 
sabin_files = [file for file in sabin_files if ".gb" in file] # Filter for files with .gb extension

# Get paths of .gb files
sabin_paths = ["../data/sabin123_reference_genomes/" + file for file in sabin_files]

# Create dictionary
sabin_files = [file.replace(".gb", "") for file in sabin_files] # Remove .gb extension
sabin_dict = dict(zip(sabin_files, sabin_paths))

# Order dict to start with Sabin1
sabin_dict = dict(sorted(sabin_dict.items()))

# Extract gene sequences for Sabin strains using previously defined function
sabin_gene_dict = {}
for sabin_path in sabin_dict.values():
    gene_dict = extract_gene_sequences(sabin_path)
    sabin = [key for key in sabin_dict.keys() if sabin_dict[key] == sabin_path][0]
    sabin_gene_dict[sabin] = gene_dict


## Run pairwise sequence alignment using biotite

import biotite.sequence as seq
import biotite.sequence.align as align

# Test Align module
matrix = align.SubstitutionMatrix.std_nucleotide_matrix() # Use standard nucleotide substitution matrix
alignment = align.align_optimal(seq.NucleotideSequence(sabin_gene_dict["Sabin1"]["VP1"]), seq.NucleotideSequence(evc_gene_dict["EnterovirusC102"]["VP1"]), matrix)[0]
pident = align.get_sequence_identity(alignment) # Calculate percentage sequence identity

# Do this for all gene sequences of all EV-C reference genomes against respective Sabin strains
# Output: dataframe with percentage sequence identity for each gene of each EV-C reference genome against respective Sabin strain
# Columns: EV-C reference genome, gene, Sabin strain, percentage sequence identity, length of gene in EV-C reference genome

# Create empty dataframe
import pandas as pd
pident_df = pd.DataFrame(columns=["evc_ref", "gene", "sabin_strain", "pident_ref", "evc_gene_length", "sabin_gene_length"])

# Iterate over EV-C reference genomes
for evc in evc_gene_dict.keys():
    # Iterate over genes
    for gene in evc_gene_dict[evc].keys():
        # Iterate over Sabin strains
        for sabin in sabin_gene_dict.keys():
            seq_evc = seq.NucleotideSequence(evc_gene_dict[evc][gene])
            seq_sabin = seq.NucleotideSequence(sabin_gene_dict[sabin][gene])
            
            # Calculate sequence identity
            alignment = align.align_optimal(seq_sabin, seq_evc, matrix)[0]
            pident = align.get_sequence_identity(alignment)
            
            # Get length of the genes in the EV-C and Sabin reference genomes, respectively
            evc_gene_length = len(evc_gene_dict[evc][gene])
            sabin_gene_length = len(sabin_gene_dict[sabin][gene])
            
            # Append to dataframe using pd.concat
            pident_df = pd.concat([pident_df, pd.DataFrame([[evc, gene, sabin, pident, evc_gene_length, sabin_gene_length]], columns=["evc_ref", "gene", "sabin_strain", "pident_ref", "evc_gene_length", "sabin_gene_length"])], ignore_index=True)
            
# Save dataframe to csv
pident_df.to_csv("../data/pident_evc_sabin.csv", index=False)