import os
import pandas as pd
import matplotlib.pyplot as plt

# Directory containing the output CSVs
output_dir = "output"
plot_dir = "plots"

# Ensure the plot directory exists
os.makedirs(plot_dir, exist_ok=True)

# Mapping of files to their respective x-axis variable
file_mapping = {
    "lda1_experiment.csv": "lda1",
    "lda2_experiment.csv": "lda2",
    "measurements.csv": "num_measurements",
    "modules_experiment.csv": "num_modules",
    "cells_experiment.csv" : "num_cells",
}

# Load and plot each file
for filename, x_var in file_mapping.items():
    file_path = os.path.join(output_dir, filename)
    
    if os.path.exists(file_path):
        # Load data
        df = pd.read_csv(file_path)

        # Create a scatter plot
        plt.figure(figsize=(8, 6))
        for gene_set_size in sorted(df["gene_set_size"].unique()):
            subset = df[df["gene_set_size"] == gene_set_size]
            plt.scatter(subset[x_var], subset["avg_best_score"], label=f"Gene Set {gene_set_size}", alpha=0.7)

        # Customize plot
        plt.xlabel(x_var.replace("_", " ").title())
        plt.ylabel("Avg Best Score")
        plt.title(f"{x_var.replace('_', ' ').title()} vs. Avg Best Score")
        plt.legend(title="Gene Set Size")
        plt.grid(False)

        # Save plot
        plot_path = os.path.join(plot_dir, f"{filename.replace('.csv', '.png')}")
        plt.savefig(plot_path, dpi=300)
        plt.close()

        print(f"Saved plot: {plot_path}")

print("All scatter plots generated and saved in 'plots/' directory.")