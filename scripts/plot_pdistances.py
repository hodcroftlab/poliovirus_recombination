

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
sns.relplot(data=pdistances, x="dist_vp1", y="dist_3dpol", kind="scatter")

# Only looking at VP1 vs 3Dpol for now
# Filter for sequences where one of the distances is much higher than the other

# Calculate the difference of the two distances and take absolute value
pdistances["diff_vp1_3dpol"] = abs(pdistances["dist_vp1"] - pdistances["dist_3dpol"])

# Sort by difference (descending)
pdistances = pdistances.sort_values("diff_vp1_3dpol", ascending=False)

