import numpy as np
import pandas as pd
import spams  # Ensure this is imported
from simulation import run_simulation
from itertools import product

THREADS = 10  # Ensure it's globally defined as an integer

# Set fixed parameters
adata_path = "dataset/pmotorcortex/pmotorcortex.h5ad"
gene_set_size = 1000
num_cells = 10000
dataset_dir = "./dataset/pmotorcortex/pmotorcortex_mouse"
num_runs = 10  # Reduce initially for testing

# Define parameter ranges
measurement_values = list(range(10, 101, 10)) + list(range(150, 501, 50))
num_modules_values = np.linspace(10, 100, 10, dtype=int)
lda1_values = np.linspace(1, 10, 10)
lda2_values = np.linspace(0.1, 1, 10)
sparsity_values = np.linspace(0.01, 0.1, 10)

# Prepare results storage
results = []

# Iterate over all parameter combinations
for num_measurements, num_modules, lda1, lda2, sparsity in product(
    measurement_values, num_modules_values, lda1_values, lda2_values, sparsity_values
):
    print(f"Running simulation for measurements={num_measurements}, modules={num_modules}, lda1={lda1}, lda2={lda2}, sparsity={sparsity}")

    # Ensure correct data types
    num_measurements = int(num_measurements)
    num_modules = int(num_modules)
    lda1 = float(lda1)
    lda2 = float(lda2)
    sparsity = float(sparsity)

    total_score = 0
    for _ in range(num_runs):
        res = run_simulation(
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
        total_score += res["best_score"]
        del res  # Free memory

    avg_score = total_score / num_runs
    results.append([num_measurements, num_modules, lda1, lda2, sparsity, avg_score])

# Convert results to DataFrame
results_df = pd.DataFrame(
    results, columns=["num_measurements", "num_modules", "lda1", "lda2", "sparsity", "avg_score"]
)

# Save results
results_df.to_csv("simulation_results.csv", index=False)

# Display results
import ace_tools as tools
tools.display_dataframe_to_user(name="Simulation Results", dataframe=results_df)