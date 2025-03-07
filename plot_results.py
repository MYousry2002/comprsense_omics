import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the final results
df = pd.read_csv("results/grid_search/final_results.csv")

# Create output directory
os.makedirs("results/grid_search/plots", exist_ok=True)

# Define parameters to analyze
parameters = ["num_measurements", "num_modules", "lda1", "lda2", "sparsity", "gene_set_size"]

# Generate plots for each parameter
for param in parameters:
    plt.figure(figsize=(10, 6))
    plt.scatter(df[param], df["avg_best_score"], alpha=0.5)
    plt.xlabel(param)
    plt.ylabel("Average Best Score")
    plt.title(f"Effect of {param} on Best Score")
    plt.grid(True)
    
    # Save the plot
    plot_path = f"results/grid_search/plots/{param}_vs_best_score.png"
    plt.savefig(plot_path)
    print(f"Plot saved: {plot_path}")