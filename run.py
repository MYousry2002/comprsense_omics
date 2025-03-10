import os
import sys
import numpy as np
import pandas as pd
import spams
from simulation import run_simulation
import itertools

THREADS = 10  # Ensure it's globally defined as an integer

# Ensure results directory exists
os.makedirs("results/grid_search", exist_ok=True)

# Set fixed parameters
adata_path = "dataset/pmotorcortex/pmotorcortex.h5ad"
num_cells = 10000
dataset_dir = "./dataset/pmotorcortex/pmotorcortex_mouse"
num_repeats = 10  # Ensure each parameterization runs 100 times  # Number of times to repeat each simulation

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

# Read job array index (allow manual override for local testing)
job_id = int(os.getenv('SGE_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '1')) - 1

# Ensure valid job ID
if job_id >= len(parameter_combinations):
    print(f"Job ID {job_id} out of range. Exiting.")
    exit(1)

# Extract parameters for this job
num_measurements, num_modules, lda1, lda2, sparsity, gene_set_size = parameter_combinations[job_id]
print(f"Running simulation with: num_measurements={num_measurements}, num_modules={num_modules}, lda1={lda1}, lda2={lda2}, sparsity={sparsity}, gene_set_size={gene_set_size}")

# Ensure correct data types
num_measurements = int(num_measurements)
num_modules = int(num_modules)
lda1 = float(lda1)
lda2 = float(lda2)
sparsity = float(sparsity)

# Run simulations and take the average best score
best_scores = []
for i in range(num_repeats):
    print(f"Running repeat {i+1}/{num_repeats}...")
    results = run_simulation(
        adata_path=adata_path,
        gene_set_size=gene_set_size,
        num_cells=num_cells,
        num_measurements=num_measurements,
        min_pools_per_gene=4,
        max_pools_per_gene=4,
        sparsity=sparsity,
        num_modules=num_modules,
        lda1=lda1,
        lda2=lda2,
        dataset_dir=dataset_dir
    )
    score = results["best_score"]
    print(f"Simulation {i+1}/{num_repeats} score: {score}")
    best_scores.append(score)
    del results  # Free memory

# Compute average best score
avg_best_score = np.mean(best_scores)

# Save results
results_df = pd.DataFrame([{
    "num_measurements": num_measurements,
    "num_modules": num_modules,
    "lda1": lda1,
    "lda2": lda2,
    "sparsity": sparsity,
    "gene_set_size": gene_set_size,
    "avg_best_score": avg_best_score
}])

csv_path = "results/grid_search/all_results.csv"
results_df.to_csv(csv_path, mode='a', header=not os.path.exists(csv_path), index=False)
print(f"Results saved: {csv_path}")
