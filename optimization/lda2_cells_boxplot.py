import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the results
csv_path = "output/lda2_cells_experiment.csv"  # Update path if needed
df = pd.read_csv(csv_path)

# Ensure the plots directory exists
plot_dir = "plots"
os.makedirs(plot_dir, exist_ok=True)

# Define custom colors for lda2 values
custom_palette = {
    "0.2": "#4c72b0",   # Muted blue
    "0.8": "#dd8452",   # Muted orange
    "2.0": "#55a868",     # Muted green
    "4.0": "#c44e52",     # Muted red
    "8.0": "#8172b3",     # Muted purple
    "12.0": "#937860",    # Muted brown
    "16.0": "#da8bc3",    # Muted pink
    "20.0": "#8c8c8c"     # Muted gray
}

# Convert lda2 to categorical for plotting
df["lda2"] = df["lda2"].astype(str)  # Ensure proper categorical grouping

# Set figure size
plt.figure(figsize=(12, 6))

# Create a box plot with num_cells on the x-axis, colored by lda2
sns.boxplot(x="num_cells", y="best_score", hue="lda2", data=df, palette=custom_palette)

# Customize plot
plt.xlabel("Number of Cells", fontsize=12)
plt.ylabel("Score", fontsize=12)
plt.title("Box Plot of Best Scores Across Number of Cells and LDA2 Values (Repeats Included)", fontsize=14)
plt.legend(title="LDA2 Value")
plt.grid(True)

# Save plot
plot_path = os.path.join(plot_dir, "boxplot_lda2_cells.png")
plt.savefig(plot_path, dpi=300)
plt.show()

print(f"Box plot saved at: {plot_path}")