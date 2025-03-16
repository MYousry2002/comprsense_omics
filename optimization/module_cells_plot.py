import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the results
csv_path = "output/modules_cells_experiment.csv"
df = pd.read_csv(csv_path)

# Ensure the plots directory exists
plot_dir = "plots"
os.makedirs(plot_dir, exist_ok=True)

# Define custom colors for num_modules
custom_palette = {
    40: "#4c72b0",   # Muted blue
    70: "#dd8452",   # Muted orange
    100: "#55a868"   # Muted green
}

# Set figure size
plt.figure(figsize=(10, 6))

# Create a box plot with num_cells on the x-axis, colored by num_modules
sns.boxplot(x="num_cells", y="best_score", hue="num_modules", data=df, palette=custom_palette)

# Customize plot
plt.xlabel("Number of Cells", fontsize=12)
plt.ylabel("Best Score", fontsize=12)
plt.title("Box Plot of Best Scores Across Number of Cells and Modules", fontsize=14)
plt.legend(title="Number of Modules")
plt.grid(True)

# Save plot
plot_path = os.path.join(plot_dir, "boxplot_cells_modules.png")
plt.savefig(plot_path, dpi=300)
plt.show()

print(f"Box plot saved at: {plot_path}")