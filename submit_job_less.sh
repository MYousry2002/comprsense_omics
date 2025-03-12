#!/bin/bash

#$ -N simulation_less
#$ -l h_rt=12:00:00       # Request 12 hours runtime
#$ -pe omp 32             # Request 18 cores per node
#$ -t 1-24            # Run 300 parallel jobs
#$ -o /projectnb/cisi/myousry/comprsense_omics/logs/simulation_less.log
#$ -e /projectnb/cisi/myousry/comprsense_omics/logs/simulation_less.err

# Load environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate comprsense

# Move to working directory
cd /projectnb/cisi/myousry/comprsense_omics/ || { echo "Error: Directory not found"; exit 1; }

# Run the Python script with job array index
python run_less.py