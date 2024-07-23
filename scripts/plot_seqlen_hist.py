
### Plot the distribution of sequence lengths in a fasta file as a histogram

import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import pandas as pd
import numpy as np


def plot_seqlen_hist(fasta_file, title="Distribution of sequence lengths", binwidth=200, xstep=400):
    
    # Create dataframe of sequence lengths
    seqlen_df = pd.DataFrame(columns=["seq_id", "seq_len"])
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqlen_df = pd.concat([seqlen_df, pd.DataFrame({"seq_id": [record.id], "seq_len": [len(record.seq)]})], ignore_index=True)
            
    # Set xlimit: round max sequence length up to the nearest integer divisible by xstep
    xlimit = max(seqlen_df["seq_len"]) + (xstep - max(seqlen_df["seq_len"]) % xstep)
    
    # Plot histogram, customize plot
    plt.figure(figsize=(10, 6))
    sns.histplot(seqlen_df, x="seq_len", binwidth=binwidth, binrange=(0, xlimit))
    plt.title(title, loc = "left")
    plt.title("Nr. of sequences: " + str(len(seqlen_df)), loc = "right")
    plt.xlabel("Sequence length")
    plt.grid(axis="both", linestyle="--", alpha=0.6)
    plt.xticks(np.arange(0, xlimit, xstep), rotation=45)

# Run function for polio dataset
plot_seqlen_hist(fasta_file=snakemake.input.fasta_polio, title="Distribution of sequence lengths in the polio dataset", binwidth=snakemake.params.binwidth, xstep=snakemake.params.xstep)
plt.savefig(snakemake.output.hist_polio) # Save output plot 

# Run function for non-polio enterovirus C dataset
plot_seqlen_hist(fasta_file=snakemake.input.fasta_npevc, title="Distribution of sequence lengths in the non-polio EV-C dataset", binwidth=snakemake.params.binwidth, xstep=snakemake.params.xstep)
plt.savefig(snakemake.output.hist_npevc) # Save output plot
