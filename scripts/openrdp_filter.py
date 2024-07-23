import pandas as pd

# Read OpenRDP results from csv file
openrdp_results = pd.read_csv("../data/openrdp_results_test.csv")

# Filter pvalues < 0.05
openrdp_results_filtered = openrdp_results[openrdp_results["Pvalue"] < 0.05]

# Filter for rows where either of the parents is one of the polioviruses, and the recombinant is not
openrdp_results_filtered = openrdp_results_filtered[openrdp_results_filtered["Parent1"].str.contains("poliovirus") != openrdp_results_filtered["Parent2"].str.contains("poliovirus")]
openrdp_results_filtered = openrdp_results_filtered[~openrdp_results_filtered["Recombinant"].str.contains("poliovirus")]

