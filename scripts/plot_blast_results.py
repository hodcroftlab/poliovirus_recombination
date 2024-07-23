
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load blast results
blast_results_polio = pd.read_csv("../data/blast_results/blast_results_polio_filtered.csv")

# Plot bitscore vs. sequence identity
plt.figure(figsize=(10, 6))
plot = sns.scatterplot(data=blast_results_polio, x="pident", y="length", hue="bitscore", palette="viridis", alpha=0.8, linewidth=0)
sns.move_legend(plot, "center right")

