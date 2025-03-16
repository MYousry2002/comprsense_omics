import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = "../results/simulation/all_results.csv"
data = pd.read_csv(file_path)

# Create a box plot for lda1 vs avg_best_score
plt.figure(figsize=(8, 6))
data.boxplot(column='avg_best_score', by='lda1', grid=False)
plt.xlabel("lda1")
plt.ylabel("Average Score")
plt.title("LDA1 vs Average Score")
plt.suptitle("")  # Remove default title
plt.xticks(rotation=45)
output_filename = "lda1_vs_avg_best_score_boxplot.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
plt.show()
print(f"Plot saved as {output_filename} in the current directory.")

# Create a box plot for lda2 vs avg_best_score
plt.figure(figsize=(8, 6))
data.boxplot(column='avg_best_score', by='lda2', grid=False)
plt.xlabel("lda2")
plt.ylabel("Average Score")
plt.title("LDA2 vs Average Score")
plt.suptitle("")  # Remove default title
plt.xticks(rotation=45)
output_filename = "lda2_vs_avg_best_score_boxplot.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
plt.show()
print(f"Plot saved as {output_filename} in the current directory.")

# Create a box plot for sparsity vs avg_best_score
plt.figure(figsize=(8, 6))
data.boxplot(column='avg_best_score', by='sparsity', grid=False)
plt.xlabel("sparsity")
plt.ylabel("avg_best_score")
plt.title("Sparsity vs Avg Best Score")
plt.suptitle("")  # Remove default title
plt.xticks(rotation=45)
output_filename = "sparsity_vs_avg_best_score_boxplot.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
plt.show()
print(f"Plot saved as {output_filename} in the current directory.")

# Create a box plot for num_modules vs avg_best_score
plt.figure(figsize=(8, 6))
data.boxplot(column='avg_best_score', by='num_modules', grid=False)
plt.xlabel("num_modules")
plt.ylabel("avg_best_score")
plt.title("Num Modules vs Avg Best Score")
plt.suptitle("")  # Remove default title
plt.xticks(rotation=45)
output_filename = "num_modules_vs_avg_best_score_boxplot.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
plt.show()
print(f"Plot saved as {output_filename} in the current directory.")

# Find the row with the highest avg_best_score
highest_score_row = data.loc[data['avg_best_score'].idxmax()]
print(highest_score_row)

# Find the top 10 rows with the highest avg_best_score along with corresponding parameter values
top_10_highest_scores = data.nlargest(10, 'avg_best_score')
print(top_10_highest_scores)