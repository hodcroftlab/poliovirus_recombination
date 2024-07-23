import pandas as pd

# Load NPEVC metadata
npevc_metadata = pd.read_csv("../data/metadata/npevc_metadata.csv")

# Filter for sequences > 300 bp
npevc_metadata_filtered = npevc_metadata[npevc_metadata["Length"] > 300]

# Count the number of sequences that contain specific keywords in their title, per serotype
def keyword_percent_inout(serotype, keywords):
    # Filter for specific serotype
    serotype_metadata = npevc_metadata_filtered[npevc_metadata_filtered["serotype"] == serotype]
    
    # Count the number of sequences that contain any of the specific keywords in their title 
    df_filtered = serotype_metadata[serotype_metadata["GenBank_Title"].str.contains("|".join(keywords), case=False)]
    accn_list = df_filtered["Accession"].tolist()
    
    # Read NP-EVC BLAST results
    npevc_blast_results = pd.read_csv("../data/blast_results/blast_results_npevc_filtered.csv")

    # How many of the sequences in accn_list appear in the BLAST results?
    npevc_blast_results_filtered = npevc_blast_results[npevc_blast_results["qaccn"].isin(accn_list)]

    # Count unique appearances of these accessions
    nunique = npevc_blast_results_filtered["qaccn"].nunique()
    
    # Calculate percentage
    percentage = nunique / len(accn_list) * 100
    
    return percentage


keyword_percent_inout(serotype="Coxsackievirus A24", keywords=["complete", "3D"])
