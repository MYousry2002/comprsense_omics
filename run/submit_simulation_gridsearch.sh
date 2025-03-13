#!/bin/bash

#$ -N grid_search
#$ -l h_rt=48:00:00       # Request 48 hours runtime
#$ -pe omp 10             # Request 10 cores per node
#$ -t 1-4500              # Run 4500 parallel jobs
#$ -o /projectnb/cisi/myousry/comprsense_omics/logs/grid_search.log
#$ -e /projectnb/cisi/myousry/comprsense_omics/logs/grid_search.err

# Load environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate comprsense

# Move to working directory
cd /projectnb/cisi/myousry/comprsense_omics/ || { echo "Error: Directory not found"; exit 1; }

# Run the Python script with job array index
python simulation_gridsearch.py