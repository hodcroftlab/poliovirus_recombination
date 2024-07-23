import pandas as pd

#------------------------------------------------------------
### Polio 
#------------------------------------------------------------

## Read and process the blast output and sequence metadata

blast_results_polio = pd.read_csv(snakemake.input.blast_results_polio, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov"])
polio_metadata = pd.read_csv(snakemake.input.polio_metadata) # Query: polio sequences
npevc_metadata = pd.read_csv(snakemake.input.npevc_metadata) # Reference: non-polio enterovirus C sequences

# Select relevant columns in the metadata tables and rename them
polio_metadata = polio_metadata[["Accession", "serotype", "Length", "Nuc_Completeness", "Country", "Host", "Isolation_Source", "Collection_Date", "GenBank_Title", "Submitters"]]
polio_metadata = polio_metadata.rename(columns={"Accession": "qaccn", "serotype": "qserotype", "Length": "qsize", "Nuc_Completeness": "qcompleteness", "Country": "qcountry", "Host": "qhost", "Isolation_Source": "qsource", "Collection_Date": "qcollectiondate", "GenBank_Title": "qgbtitle", "Submitters": "qsubmitters"})

npevc_metadata = npevc_metadata[["Accession", "serotype", "Length", "Nuc_Completeness", "Country", "Host", "Isolation_Source", "Collection_Date", "GenBank_Title", "Submitters"]]
npevc_metadata = npevc_metadata.rename(columns={"Accession": "saccn", "serotype": "sserotype", "Length": "ssize", "Nuc_Completeness": "scompleteness", "Country": "scountry", "Host": "shost", "Isolation_Source": "ssource", "Collection_Date": "scollectiondate", "GenBank_Title": "sgbtitle", "Submitters": "ssubmitters"})

# Rename qseqid and sseqid columns to match metadata
blast_results_polio = blast_results_polio.rename(columns={"qseqid": "qaccn", "sseqid": "saccn"})

# Inner merge the three dataframes --> this will remove rows that have been filtered out in the metadata dataframes previously
blast_results_polio = blast_results_polio.merge(polio_metadata, on="qaccn", how="inner")
blast_results_polio = blast_results_polio.merge(npevc_metadata, on = "saccn", how="inner")


## Filter the results

min_match_length = snakemake.params.min_match_length
min_sequence_identity = snakemake.params.min_sequence_identity

# Filter by match length (try at least 300)
blast_results_polio = blast_results_polio[blast_results_polio["length"] >= min_match_length]

# Filter by sequence identity (try 85)
blast_results_polio = blast_results_polio[blast_results_polio["pident"] >= min_sequence_identity]

# Sort by bitscore
blast_results_polio = blast_results_polio.sort_values(by="bitscore", ascending=False)

# Order columns to have query and sequence genome names, accession numbers and serotypes first
blast_results_polio = blast_results_polio[["qaccn", "qserotype", "qgbtitle", "saccn", "sserotype", "sgbtitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov", "qsize", "qcompleteness", "qcountry", "qhost", "qsource", "qcollectiondate", "qsubmitters", "ssize", "scompleteness", "scountry", "shost", "ssource", "scollectiondate", "ssubmitters"]]

# The majority of very high sequence identity matches are with unidentified enteroviruses C
# These are likely to be polio
# Filter out these matches
blast_results_polio = blast_results_polio[~(blast_results_polio["sserotype"] == "Unidentified Enterovirus C")]

# Save filtered dataframe as .csv
blast_results_polio.to_csv(snakemake.output.filtered_blast_polio, index=False)
blast_results_polio_unique_accn = blast_results_polio.drop_duplicates(subset="qaccn").to_csv(snakemake.output.unique_accn_blast_polio, index=False) # Filter for first appearance of a query sequence
blast_results_polio_unique_submitters = blast_results_polio.drop_duplicates(subset="qsubmitters").to_csv(snakemake.output.unique_submitter_blast_polio, index=False) # Filter for first appearance of a certain (set of) submitter(s)

#------------------------------------------------------------
### Non-polio enterovirus C
#------------------------------------------------------------

## Read and process the blast output and sequence metadata
blast_results_npevc = pd.read_csv(snakemake.input.blast_results_npevc, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov"])

# Rename columns in the NPEVC metadata: this is the query metadata now
npevc_metadata = npevc_metadata.rename(columns={"saccn": "qaccn", "sserotype": "qserotype", "ssize": "qsize", "scompleteness": "qcompleteness", "scountry": "qcountry", "shost": "qhost", "ssource": "qsource", "scollectiondate": "qcollectiondate", "sgbtitle": "qgbtitle", "ssubmitters": "qsubmitters"})

# Read sabin123_metadata.csv
sabin123_metadata = pd.read_csv(snakemake.input.sabin123_metadata)
sabin123_metadata = sabin123_metadata[["Accession", "Organism_Name", "Length", "Release_Date", "Nuc_Completeness", "GenBank_Title"]]
sabin123_metadata["Organism_Name"] = sabin123_metadata["Organism_Name"].str.replace("Human poliovirus", "Sabin") # Map serotype from organism name
sabin123_metadata = sabin123_metadata.rename(columns={"Accession": "saccn", "Organism_Name": "sserotype", "Length": "ssize", "Release_Date": "sreleasedate", "Nuc_Completeness": "scompleteness", "GenBank_Title": "sgbtitle"})

# Rename qseqid and sseqid columns to match metadata
blast_results_npevc = blast_results_npevc.rename(columns={"qseqid": "qaccn", "sseqid": "saccn"})

# Inner merge the three dataframes --> this will remove rows that have been filtered out in the metadata dataframes previously
blast_results_npevc = blast_results_npevc.merge(npevc_metadata, on="qaccn", how="inner")
blast_results_npevc = blast_results_npevc.merge(sabin123_metadata, on="saccn", how="inner")


## Filter the results by match length and sequence identity
blast_results_npevc = blast_results_npevc[blast_results_npevc["length"] >= min_match_length]
blast_results_npevc = blast_results_npevc[blast_results_npevc["pident"] >= min_sequence_identity]

# Sort by bitscore
blast_results_npevc = blast_results_npevc.sort_values(by="bitscore", ascending=False)

# Order columns to have query and sequence genome names, accession numbers and serotypes first
blast_results_npevc = blast_results_npevc[["qaccn", "qserotype", "qgbtitle", "saccn", "sserotype", "sgbtitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov", "qsize", "qcompleteness", "qcountry", "qhost", "qsource", "qcollectiondate", "qsubmitters", "ssize", "sreleasedate", "scompleteness"]]

# Save filtered dataframe as .csv
blast_results_npevc.to_csv(snakemake.output.filtered_blast_npevc, index=False)
blast_results_npevc_unique_accn = blast_results_npevc.drop_duplicates(subset="qaccn").to_csv(snakemake.output.unique_accn_blast_npevc, index=False) # Filter for first appearance of a query sequence
blast_results_npevc_unique_submitters = blast_results_npevc.drop_duplicates(subset="qsubmitters").to_csv(snakemake.output.unique_submitter_blast_npevc, index=False) # Filter for first appearance of a certain (set of) submitter(s)
