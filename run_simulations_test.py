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


# Set parameters
adata_path = "dataset/pmotorcortex/pmotorcortex.h5ad"  # Path to your AnnData file
gene_set_size = 1000  # Choose from 500, 1000, or 5000
num_cells = 10000  # Number of cells to include in the simulation
num_measurements = 200  # Number of measurement pools
min_pools_per_gene = 4  # Minimum number of pools per gene
max_pools_per_gene = 4  # Maximum number of pools per gene
sparsity = 0.02  # Sparsity constraint for sparse decoding
num_modules = 50  # Number of gene modules
lda1=8
lda2=0.2
dataset_dir = "./dataset/pmotorcortex/pmotorcortex_mouse"  # Directory where datasets are stored

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