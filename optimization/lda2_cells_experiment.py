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

# Set fixed parameters
num_measurements = 100  # Fixed at 100
gene_set_size = 1000  # Fixed at 1000
lda1 = 8.0  # Keeping lda1 fixed
sparsity = 0.02
num_repeats = 10  # Ensure each parameterization runs multiple times

# Define parameter search space (Varying lda2 and num_cells)
lda2_values = [0.2, 0.8, 2, 4, 8, 12, 16, 20]  # 8 values
num_cells_values = [1000, 10000, 20000, 50000, 75000, 100000]  # 6 values

# Generate all parameter combinations
parameter_combinations = [(lda2, num_cells) 
                          for lda2 in lda2_values
                          for num_cells in num_cells_values]

# Read job array index (allow manual override for local testing)
job_id = int(os.getenv('SGE_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '1')) - 1

# Ensure valid job ID
if job_id >= len(parameter_combinations):
    print(f"Job ID {job_id} out of range. Exiting.")
    exit(1)

# Extract parameters for this job
lda2, num_cells = parameter_combinations[job_id]
print(f"Running simulation with: num_cells={num_cells}, lda2={lda2}, num_measurements={num_measurements}, gene_set_size={gene_set_size}, lda1={lda1}, sparsity={sparsity}")

# Ensure correct data types
num_cells = int(num_cells)
lda2 = float(lda2)

# Prepare DataFrame for tracking results across repeats
results_list = []

# Run simulations and track every repeat
for repeat_id in range(1, num_repeats + 1):
    print(f"Running repeat {repeat_id}/{num_repeats}...")
    results = run_simulation(
        adata_path=os.path.join(ROOT_DIR, "dataset/pmotorcortex/pmotorcortex.h5ad"),
        gene_set_size=gene_set_size,  # Fixed at 1000
        num_cells=num_cells,  # Now varying num_cells
        num_measurements=num_measurements,  # Fixed at 100
        min_pools_per_gene=4,
        max_pools_per_gene=4,
        sparsity=sparsity,
        num_modules=40,  # Fixed at 40
        lda1=lda1,  # Fixed at 8.0
        lda2=lda2,  # Now varying lda2
        dataset_dir=os.path.join(ROOT_DIR, "dataset/pmotorcortex/pmotorcortex_mouse")
    )
    score = results["best_score"]
    print(f"Simulation {repeat_id}/{num_repeats} score: {score}")
    
    # Store each repeat separately
    results_list.append({
        "repeat_id": repeat_id,
        "num_cells": num_cells,
        "lda2": lda2,
        "num_measurements": num_measurements,  # Fixed at 100
        "gene_set_size": gene_set_size,  # Fixed at 1000
        "lda1": lda1,  # Fixed at 8.0
        "sparsity": sparsity,
        "best_score": score  # Store score per repeat
    })
    
    del results  # Free memory

# Convert results to DataFrame and save
results_df = pd.DataFrame(results_list)

csv_path = "output/lda2_cells_experiment.csv"
results_df.to_csv(csv_path, mode='a', header=not os.path.exists(csv_path), index=False)

print(f"Results saved: {csv_path}")