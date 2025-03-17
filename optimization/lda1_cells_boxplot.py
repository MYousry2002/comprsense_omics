import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the results
csv_path = "output/lda1_cells_experiment.csv"  # Update path if needed
df = pd.read_csv(csv_path)

# Ensure the plots directory exists
plot_dir = "plots"
os.makedirs(plot_dir, exist_ok=True)

# Define custom colors for lda1 values
custom_palette = {
    "2.0": "#4c72b0",   # Muted blue
    "4.0": "#dd8452",   # Muted orange
    "8.0": "#55a868",   # Muted green
    "12.0": "#c44e52",  # Muted red
    "16.0": "#8172b3",  # Muted purple
    "20.0": "#937860"   # Muted brown
}

# Convert lda1 to categorical for plotting
df["lda1"] = df["lda1"].astype(str)  # Ensure proper categorical grouping

# Set figure size
plt.figure(figsize=(12, 6))

# Create a box plot with num_cells on the x-axis, colored by lda1
sns.boxplot(x="num_cells", y="best_score", hue="lda1", data=df, palette=custom_palette)

# Customize plot
plt.xlabel("Number of Cells", fontsize=12)
plt.ylabel("Score", fontsize=12)
plt.title("Box Plot of Best Scores Across Number of Cells and LDA1 Values (Repeats Included)", fontsize=14)
plt.legend(title="LDA1 Value")
plt.grid(True)

# Save plot
plot_path = os.path.join(plot_dir, "boxplot_lda1_cells.png")
plt.savefig(plot_path, dpi=300)
plt.show()

print(f"Box plot saved at: {plot_path}")