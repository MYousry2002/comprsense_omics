#!/bin/bash

#$ -N comprsense_simulation
#$ -l h_rt=24:00:00
#$ -pe omp 10
#$ -o /projectnb/cisi/myousry/comprsense_omics/comprsense_simulation.log
#$ -e /projectnb/cisi/myousry/comprsense_omics/comprsense_simulation.err

# Ensure Conda is sourced correctly from your installation path
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate comprsense

# Move to the working directory
cd /projectnb/cisi/myousry/comprsense_omics/ || { echo "Error: Directory not found"; exit 1; }

# Check if the script exists before running
if [[ ! -f run_simulations.py ]]; then
    echo "Error: run_simulations.py not found!"
    exit 1
fi

# Run the simulation script
python run_simulations.py


#!/bin/bash

