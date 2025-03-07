import pandas as pd
import glob

# Load all result files
csv_files = glob.glob("results/grid_search/results_*.csv")
df_list = [pd.read_csv(f) for f in csv_files]

# Merge results
df_results = pd.concat(df_list, ignore_index=True)

# Save final combined results
df_results.to_csv("results/grid_search/final_results.csv", index=False)
print("Final results saved as: results/grid_search/final_results.csv")