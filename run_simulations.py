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

THREADS = 10

from simulation import run_simulation


# Define a list of values for num_measurements to search over
num_measurements_values = list(range(10, 100, 10)) + list(range(100, 500, 50)) + [500]

# Placeholder for results storage
simulation_results = []

# Run the simulation for each num_measurements value
for num_measurements in num_measurements_values:
    results = run_simulation(
        adata_path = "dataset/cerebellum/cb_adult_mouse.h5ad",
        gene_set_size=1000,
        num_cells=10000,
        num_measurements=num_measurements,
        min_pools_per_gene=10,
        max_pools_per_gene=10,
        sparsity=0.02,
        num_modules=50,
        dataset_dir = "./dataset/cerebellum/cb_mouse"
    )
    
    # Store results
    simulation_results.append({
        "num_measurements": num_measurements,
        "best_matrix_shape": results["best_matrix"].shape,
        "best_score": results["best_score"],
        "best_coherence_scores": results["best_coherence_scores"][:5],
        "best_reconstruction_scores": results["best_reconstruction_scores"][:5]
    })

# Convert results to a DataFrame and save or display
df_results = pd.DataFrame(simulation_results)
print(df_results)

# save to a CSV file
df_results.to_csv("results/measurements_simulation/measurements_simulation_results.csv", index=False)