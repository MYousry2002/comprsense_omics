import os
import sys
import numpy as np
import pandas as pd
import spams
from simulation import run_simulation

THREADS = 10  # Ensure it's globally defined as an integer

# Ensure results directory exists
os.makedirs("results/simulation", exist_ok=True)

# Set fixed parameters
adata_path = "dataset/pmotorcortex/pmotorcortex.h5ad"
num_cells = 10000
dataset_dir = "./dataset/pmotorcortex/pmotorcortex_mouse"
num_modules = 40
lda1 = 8.0
lda2 = 0.8
sparsity = 0.02
num_repeats = 10  # Ensure each parameterization runs multiple times

# Define parameter search space with remaining variable parameters
num_measurements_values = [10, 25, 50, 75, 100, 200, 350, 500]  # 8 values
gene_set_sizes = [500, 1000, 5000]  # 3 values

# Generate all parameter combinations
parameter_combinations = [(num_measurements, gene_set_size) 
                          for num_measurements in num_measurements_values 
                          for gene_set_size in gene_set_sizes]

# Read job array index (allow manual override for local testing)
job_id = int(os.getenv('SGE_TASK_ID', sys.argv[1] if len(sys.argv) > 1 else '1')) - 1

# Ensure valid job ID
if job_id >= len(parameter_combinations):
    print(f"Job ID {job_id} out of range. Exiting.")
    exit(1)

# Extract parameters for this job
num_measurements, gene_set_size = parameter_combinations[job_id]
print(f"Running simulation with: num_measurements={num_measurements}, num_modules={num_modules}, lda1={lda1}, lda2={lda2}, sparsity={sparsity}, gene_set_size={gene_set_size}")

# Ensure correct data types
num_measurements = int(num_measurements)
gene_set_size = int(gene_set_size)

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

csv_path = "results/simulation/less_results.csv"
results_df.to_csv(csv_path, mode='a', header=not os.path.exists(csv_path), index=False)
print(f"Results saved: {csv_path}")