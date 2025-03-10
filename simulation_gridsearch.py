import os
import numpy as np
import pandas as pd
import itertools
from simulation import run_simulation

# Define parameter search space
num_measurements_values = list(range(10, 100, 10)) + list(range(100, 500, 50)) + [500]
num_modules_values = np.linspace(10, 100, 10, dtype=int)
lda1_values = np.linspace(1, 20, 10)
lda2_values = np.linspace(0.1, 2, 10)
sparsity_values = np.linspace(0.01, 0.1, 10)
gene_set_sizes = [500, 1000, 5000]

# Generate all parameter combinations
parameter_combinations = list(itertools.product(
    num_measurements_values, num_modules_values, lda1_values, lda2_values, sparsity_values, gene_set_sizes
))

# Read job array index
job_id = int(os.getenv('SGE_TASK_ID', '1')) - 1

# Ensure valid job ID
if job_id >= len(parameter_combinations):
    print(f"Job ID {job_id} out of range. Exiting.")
    exit(1)

# Extract parameters for this job
num_measurements, num_modules, lda1, lda2, sparsity, gene_set_size = parameter_combinations[job_id]
print(f"Running simulation with: num_measurements={num_measurements}, num_modules={num_modules}, lda1={lda1}, lda2={lda2}, sparsity={sparsity}, gene_set_size={gene_set_size}")

# Run ~100 simulations and take the average best score
num_repeats = 100
best_scores = []

for _ in range(num_repeats):
    results = run_simulation(
        adata_path="dataset/pmotorcortex/pmotorcortex.h5ad",
        gene_set_size=gene_set_size,
        num_cells=10000,
        num_measurements=num_measurements,
        min_pools_per_gene=4,
        max_pools_per_gene=4,
        sparsity=sparsity,
        num_modules=num_modules,
        lda1=lda1,
        lda2=lda2,
        dataset_dir="./dataset/pmotorcortex/pmotorcortex_mouse"
    )
    best_scores.append(results["best_score"])

# Compute average best score
avg_best_score = np.mean(best_scores)

# Save results
df_results = pd.DataFrame([{
    "num_measurements": num_measurements,
    "num_modules": num_modules,
    "lda1": lda1,
    "lda2": lda2,
    "sparsity": sparsity,
    "gene_set_size": gene_set_size,
    "avg_best_score": avg_best_score
}])

csv_path = f"results/grid_search/results_{job_id}.csv"
df_results.to_csv(csv_path, index=False)
print(f"Results saved: {csv_path}")