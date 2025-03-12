import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = "../results/simulation/less_results_test.csv"
data = pd.read_csv(file_path)

# Create a scatter plot for lda1 vs avg_best_score
plt.figure(figsize=(8, 6))
plt.scatter(data['num_measurements'], data['avg_best_score'], alpha=0.7)
plt.xlabel("num_measurements")
plt.ylabel("avg_best_score")
plt.title("num_measurements vs Avg Best Score")
plt.grid(False)

# Save the plot
output_filename = "num_measurements_vs_avg_best_score.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Plot saved as {output_filename} in the current directory.")


# Define unique gene set sizes
gene_sizes = [500, 1000, 5000]

# Map gene_set_size to colors (using categorical colormap)
color_map = {500: "blue", 1000: "green", 5000: "red"}  
colors = data["gene_set_size"].map(color_map)

# Create scatter plot
plt.figure(figsize=(8, 6))
scatter = plt.scatter(data['num_measurements'], data['avg_best_score'], c=colors, alpha=0.7)

# Add legend manually
legend_labels = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[g], markersize=8, label=str(g)) for g in gene_sizes]
plt.legend(handles=legend_labels, title="Gene Set Size")

# Labels and title
plt.xlabel("num_measurements")
plt.ylabel("avg_best_score")
plt.title("num_measurements vs Avg Best Score (Colored by Gene Set Size)")
plt.grid(False)

# Save the plot
output_filename = "num_measurements_vs_avg_best_score_genesetsize_colored.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Plot saved as {output_filename} in the current directory.")