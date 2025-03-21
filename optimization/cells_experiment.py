#!/usr/bin/env python

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

# Fixed num_measurements, now varying num_cells
num_measurements = 200  # This is now fixed
num_modules = 40
lda1 = 8.0
lda2 = 10.0
sparsity = 0.02
num_repeats = 10  # Ensure each parameterization runs multiple times

# Define parameter search space (Varying num_cells instead of num_measurements)
num_cells_values = [1000, 5000, 10000, 20000, 30000, 50000, 75000, 100000]  # 8 values
gene_set_sizes = [500, 1000, 5000]  # 3 values

# Generate all parameter combinations
parameter_combinations = [(num_cells, gene_set_size) 
                          for num_cells in num_cells_values 
                          for gene_set_size in gene_set_sizes]

# Read job array index (allow manual override for local testing)
job_id = int(os.getenv('SGE_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '1')) - 1

# Ensure valid job ID
if job_id >= len(parameter_combinations):
    print(f"Job ID {job_id} out of range. Exiting.")
    exit(1)

# Extract parameters for this job
num_cells, gene_set_size = parameter_combinations[job_id]
print(f"Running simulation with: num_cells={num_cells}, num_measurements={num_measurements}, num_modules={num_modules}, lda1={lda1}, lda2={lda2}, sparsity={sparsity}, gene_set_size={gene_set_size}")

# Ensure correct data types
num_cells = int(num_cells)
gene_set_size = int(gene_set_size)

# Prepare a list to store repeat results
results_list = []

# Run simulations and record all repeats
for repeat_id in range(1, num_repeats + 1):
    print(f"Running repeat {repeat_id}/{num_repeats}...")
    results = run_simulation(
        adata_path=adata_path,
        gene_set_size=gene_set_size,
        num_cells=num_cells,  # Now varying num_cells
        num_measurements=num_measurements,  # Fixed at 100
        min_pools_per_gene=4,
        max_pools_per_gene=4,
        sparsity=sparsity,
        num_modules=num_modules,
        lda1=lda1,
        lda2=lda2,
        dataset_dir=dataset_dir
    )
    
    score = results["best_score"]
    print(f"Simulation {repeat_id}/{num_repeats} score: {score}")

    # Store individual repeat results
    results_list.append({
        "repeat_id": repeat_id,
        "num_cells": num_cells,
        "num_modules": num_modules,
        "num_measurements": num_measurements,
        "gene_set_size": gene_set_size,
        "lda1": lda1,
        "lda2": lda2,
        "sparsity": sparsity,
        "best_score": score  # Store each score separately
    })

    del results  # Free memory

# Convert to DataFrame and save each repeat separately
results_df = pd.DataFrame(results_list)

csv_path = "output/cells_experiment_repeats.csv"
results_df.to_csv(csv_path, mode='a', header=not os.path.exists(csv_path), index=False)

print(f"Results saved: {csv_path}")