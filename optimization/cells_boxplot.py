import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the results
csv_path = "output/cells_experiment_repeats.csv"  # Update path if needed
df = pd.read_csv(csv_path)

# Ensure the plots directory exists
plot_dir = "plots"
os.makedirs(plot_dir, exist_ok=True)

# Define custom colors for gene_set_size
custom_palette = {
    500: "#4c72b0",   # Muted blue
    1000: "#dd8452",  # Muted orange
    5000: "#55a868"   # Muted green
}

# Set figure size
plt.figure(figsize=(12, 6))

# Create a box plot with num_cells on the x-axis, colored by gene_set_size
sns.boxplot(x="num_cells", y="best_score", hue="gene_set_size", data=df, palette=custom_palette)

# Customize plot
plt.xlabel("Number of Cells", fontsize=12)
plt.ylabel("Score", fontsize=12)
plt.title("Box Plot of Best Scores Across Number of Cells and Gene Set Sizes (Repeats Included)", fontsize=14)
plt.legend(title="Gene Set Size")
plt.grid(True)

# Save plot
plot_path = os.path.join(plot_dir, "boxplot_cells_gene_sets.png")
plt.savefig(plot_path, dpi=300)
plt.show()

print(f"Box plot saved at: {plot_path}")