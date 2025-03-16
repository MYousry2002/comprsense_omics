import os
import sys
import numpy as np
import pandas as pd
import spams

# Get the absolute path of the project's root directory
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Add `src/` to sys.path
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))

# Import the simulation function
from simulation import run_simulation

THREADS = 10  # Ensure it's globally defined as an integer

# Ensure results directory exists
os.makedirs(os.path.join(ROOT_DIR, "results"), exist_ok=True)

# Set fixed parameters with absolute paths
adata_path = os.path.join(ROOT_DIR, "dataset/pmotorcortex/pmotorcortex.h5ad")
dataset_dir = os.path.join(ROOT_DIR, "dataset/pmotorcortex/pmotorcortex_mouse")


# Set parameters
gene_set_size = 1000  # Choose from 500, 1000, or 5000
num_cells = 10000  # Number of cells to include in the simulation
num_measurements = 200  # Number of measurement pools
min_pools_per_gene = 4  # Minimum number of pools per gene
max_pools_per_gene = 4  # Maximum number of pools per gene
sparsity = 0.02  # Sparsity constraint for sparse decoding
num_modules = 50  # Number of gene modules
lda1=8
lda2=8

# Run the simulation
results = run_simulation(
    adata_path=adata_path,
    gene_set_size=gene_set_size,
    num_cells=num_cells,
    num_measurements=num_measurements,
    min_pools_per_gene=min_pools_per_gene,
    max_pools_per_gene=max_pools_per_gene,
    sparsity=sparsity,
    num_modules=num_modules,
    lda1=lda1,
    lda2=lda2,
    dataset_dir=dataset_dir
)

# Print results
print("Best Measurement Matrix Shape:", results["best_matrix"].shape)
print("Best Score:", results["best_score"])
print("Top 5 Coherence Scores:", results["best_coherence_scores"][:5])
print("Top 5 Reconstruction Scores:", results["best_reconstruction_scores"][:5])