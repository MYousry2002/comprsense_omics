import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = "../results/simulation/all_results.csv"
data = pd.read_csv(file_path)

# Create a scatter plot for lda1 vs avg_best_score
plt.figure(figsize=(8, 6))
plt.scatter(data['lda1'], data['avg_best_score'], alpha=0.7)
plt.xlabel("lda1")
plt.ylabel("avg_best_score")
plt.title("LDA1 vs Avg Best Score")
plt.grid(True)

# Save the plot
output_filename = "lda1_vs_avg_best_score.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Plot saved as {output_filename} in the current directory.")


# Create a scatter plot for lda1 vs avg_best_score
plt.figure(figsize=(8, 6))
plt.scatter(data['lda2'], data['avg_best_score'], alpha=0.7)
plt.xlabel("lda2")
plt.ylabel("avg_best_score")
plt.title("LDA2 vs Avg Best Score")
plt.grid(True)

# Save the plot
output_filename = "lda2_vs_avg_best_score.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()
print(f"Plot saved as {output_filename} in the current directory.")



# Create a scatter plot for lda1 vs avg_best_score
plt.figure(figsize=(8, 6))
plt.scatter(data['sparsity'], data['avg_best_score'], alpha=0.7)
plt.xlabel("sparsity")
plt.ylabel("avg_best_score")
plt.title("sparsity vs Avg Best Score")
plt.grid(True)

# Save the plot
output_filename = "sparsity_vs_avg_best_score.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Plot saved as {output_filename} in the current directory.")



# Create a scatter plot for lda1 vs avg_best_score
plt.figure(figsize=(8, 6))
plt.scatter(data['num_modules'], data['avg_best_score'], alpha=0.7)
plt.xlabel("num_modules")
plt.ylabel("avg_best_score")
plt.title("num_modules vs Avg Best Score")
plt.grid(True)

# Save the plot
output_filename = "num_modules_vs_avg_best_score.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

print(f"Plot saved as {output_filename} in the current directory.")



# Find the row with the highest avg_best_score
highest_score_row = data.loc[data['avg_best_score'].idxmax()]

# Print the result
print(highest_score_row)


# Find the top 10 rows with the highest avg_best_score along with corresponding parameter values
top_10_highest_scores = data.nlargest(10, 'avg_best_score')

# Print the result
print(top_10_highest_scores)