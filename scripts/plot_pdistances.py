

# Plot pdistance of VP1 vs. 3Dpol and VP1 vs. 5'UTR

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load the data
pdistances = pd.read_csv("../data/pdistances/pdistances.csv")

# Plot the data
sns.set_style("whitegrid")
plt.figure(figsize=(8, 6))
plt.xlabel("p-distance VP1")
plt.ylabel("p-distance 3Dpol")
plt.title("p-distance VP1 vs. 3Dpol")
plt.scatter(pdistances["pdist_VP1"], pdistances["pdist_3Dpol"], alpha=0.3)

