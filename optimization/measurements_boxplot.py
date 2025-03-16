import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the results
csv_path = "output/measurements_boxplot_experiment.csv"
df = pd.read_csv(csv_path)

# Ensure the plots directory exists
plot_dir = "plots"
os.makedirs(plot_dir, exist_ok=True)

# Define custom color palette
custom_palette = {
    500: "#4c72b0",   # Muted blue
    1000: "#dd8452",  # Muted orange
    5000: "#55a868"   # Muted green
}

# Set figure size
plt.figure(figsize=(10, 6))

# Create a box plot with the specified colors
sns.boxplot(x="num_measurements", y="best_score", hue="gene_set_size", data=df, palette=custom_palette)

# Customize plot
plt.xlabel("Number of Measurements", fontsize=12)
plt.ylabel("Score", fontsize=12)
plt.title("Box Plot of Scores Across Measurement Values and Gene Set Size", fontsize=14)
plt.legend(title="Gene Set Size")
plt.grid(True)

# Save plot
plot_path = os.path.join(plot_dir, "boxplot_measurements.png")
plt.savefig(plot_path, dpi=300)
plt.show()

print(f"Box plot saved at: {plot_path}")