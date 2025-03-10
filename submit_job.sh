#!/bin/bash

#$ -N simulation
#$ -l h_rt=96:00:00       # Request 48 hours runtime
#$ -pe omp 18             # Request 18 cores per node
#$ -t 1-45000            # Run 4500 parallel jobs
#$ -o /projectnb/cisi/myousry/comprsense_omics/logs/simulation.log
#$ -e /projectnb/cisi/myousry/comprsense_omics/logs/simulation.err

# Load environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate comprsense

# Move to working directory
cd /projectnb/cisi/myousry/comprsense_omics/ || { echo "Error: Directory not found"; exit 1; }

# Run the Python script with job array index
python run.py