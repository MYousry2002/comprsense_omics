import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import spams
from scipy.stats import entropy
from scipy.spatial import distance
from sklearn.model_selection import train_test_split
import scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import spearmanr, pearsonr, entropy
from scipy.spatial import distance
import os
import matplotlib.pyplot as plt

THREADS = 10

from simulation import run_simulation

# Ensure results directory exists
os.makedirs("results/measurements_simulation", exist_ok=True)

# Define a list of values for num_measurements to search over
num_measurements_values = list(range(10, 100, 10)) + list(range(100, 500, 50)) + [500]

# Placeholder for results storage
simulation_results = []

# Run the simulation for each num_measurements value
for num_measurements in num_measurements_values:
    results = run_simulation(
        adata_path="dataset/cerebellum/cb_adult_mouse.h5ad",
        gene_set_size=1000,
        num_cells=10000,
        num_measurements=num_measurements,
        min_pools_per_gene=4,
        max_pools_per_gene=4,
        sparsity=0.02,
        num_modules=50,
        dataset_dir="./dataset/cerebellum/cb_mouse"
    )
    
    # Store results
    simulation_results.append({
        "num_measurements": num_measurements,
        "best_score": results["best_score"]
    })

# Convert results to a DataFrame and save
df_results = pd.DataFrame(simulation_results)
print(df_results)

# Save results to CSV
csv_path = "results/measurements_simulation/best_score_results.csv"
df_results.to_csv(csv_path, index=False)
print(f"Best score results saved to {csv_path}")

# ---- Generate and Save Best Score Plot ----
plt.figure(figsize=(10, 6))
plt.plot(df_results["num_measurements"], df_results["best_score"], marker="o", linestyle="-")
plt.xlabel("Number of Measurements")
plt.ylabel("Best Score")
plt.title("Effect of Number of Measurements on Best Score")
plt.grid(True)

# Save the plot
score_plot_path = "results/measurements_simulation/best_score_vs_measurements.png"
plt.savefig(score_plot_path)
print(f"Plot saved: {score_plot_path}")