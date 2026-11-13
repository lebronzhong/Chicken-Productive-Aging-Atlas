#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Trajectory clustering.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import pandas as pd
import numpy as np
import os
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# Analysis parameters.
random_seed = 42
n_clusters = 10
n_init_times = 100
# Input file path.
file_path = "/path/to/chicken-aging/ont/shell/14.5.LOESS_Cluster_ALL/2.Cluster/All_Tissues_Cluster_cleaned_merged.txt"
# Output directory.
output_dir = "Cluster_Results_V2"
os.makedirs(output_dir, exist_ok=True)
# Read input data.
df = pd.read_csv(file_path, sep="\t", index_col=0)
# Cluster features.
kmeans = KMeans(n_clusters=n_clusters, random_state=random_seed, n_init=n_init_times)
cluster_labels = kmeans.fit_predict(df)
# Add cluster assignments.
df_clustered = df.copy()
df_clustered["cluster_raw"] = cluster_labels
# Calculate means for the original cluster assignments.
cluster_mean = df_clustered.groupby("cluster_raw").mean()
# Order clusters by mean expression and renumber from one.
cluster_order = cluster_mean.mean(axis=1).sort_values(ascending=False).index.tolist()
cluster_mapping = {
    old: new + 1 for new, old in enumerate(cluster_order)
}  # Number clusters starting from one.
df_clustered["cluster"] = df_clustered["cluster_raw"].map(cluster_mapping)
df_clustered.drop(columns="cluster_raw", inplace=True)
# Save output.
# Export all transcripts with cluster assignments.
df_clustered.to_csv(f"{output_dir}/Clustered_transcripts_with_labels.txt", sep="\t")
# Export transcript IDs for each cluster.
for cl in sorted(df_clustered["cluster"].unique()):
    genes = df_clustered[df_clustered["cluster"] == cl].index
    genes.to_series().to_csv(f"{output_dir}/Cluster_{cl}_genes.txt", index=False, header=False)
# Export cluster means and member counts.
cluster_mean_final = df_clustered.groupby("cluster").mean()
cluster_sizes = df_clustered["cluster"].value_counts().sort_index()
cluster_mean_final["Count"] = cluster_sizes
cluster_mean_final.to_csv(f"{output_dir}/Cluster_Mean_Expression_with_Counts.txt", sep="\t")
# Plot the overall cluster trajectories.
plt.figure(figsize=(10, 6))
colors = sns.color_palette("tab10", n_clusters)
for cl in range(1, n_clusters + 1):
    mean_expr = cluster_mean_final.loc[cl].drop("Count")
    plt.plot(
        mean_expr.index,
        mean_expr.values,
        label=f"Cluster {cl} (n={int(cluster_mean_final.loc[cl]['Count'])})",
        color=colors[(cl - 1) % 10],
    )  # Keep palette indices within bounds.
plt.xlabel("Time Point")
plt.ylabel("Z-score Expression")
plt.title("Cluster Expression Trends (Mean)")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig(f"{output_dir}/Cluster_Mean_Trend_Plot.png", dpi=300)
plt.close()
