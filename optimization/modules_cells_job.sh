#!/bin/bash

#$ -N modules_cells_simulation
#$ -l h_rt=4:00:00       # Request 12 hours runtime
#$ -pe omp 32             # Request 32 cores per node
#$ -t 1-200            # Run 24 parallel jobs
#$ -o /projectnb/cisi/myousry/comprsense_omics/logs/optimization/modules_cells_simulation.log
#$ -e /projectnb/cisi/myousry/comprsense_omics/logs/optimization/modules_cells_simulation.err

# Load environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate comprsense

# Move to working directory
cd /projectnb/cisi/myousry/comprsense_omics/optimization || { echo "Error: Directory not found"; exit 1; }

# Run the Python script with job array index
python modules_cells_experiment.py