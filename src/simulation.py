import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import spams
from scipy.stats import entropy
from scipy.spatial import distance
from sklearn.model_selection import train_test_split
import scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import spearmanr, pearsonr, entropy
from scipy.spatial import distance
import os

THREADS = 10

from prepare_data import downsample_adata, split_and_save_data, load_saved_data
from measurements_matrix import find_best_coherence_matrices, find_best_reconstruction_matrices
from module_dictionary import smaf


def run_simulation(
    adata_path, 
    gene_set_size,
    num_cells, 
    num_measurements, 
    min_pools_per_gene, 
    max_pools_per_gene, 
    sparsity, 
    num_modules,
    lda1=8,
    lda2=0.2,
    dataset_dir="./dataset/"
):
    """
    Runs the full simulation pipeline to optimize measurement matrices for gene expression analysis.
    
    Parameters:
        adata_path (str): Path to the AnnData file.
        gene_set_size (int): Number of genes to include in the simulation (500, 1000, or 5000).
        num_cells (int): Number of cells to include in the simulation.
        num_measurements (int): Number of measurement pools.
        min_pools_per_gene (int): Minimum pools a gene must be assigned to.
        max_pools_per_gene (int): Maximum pools a gene can be assigned to.
        sparsity (float): Sparsity constraint for sparse decoding.
        num_modules (int): Number of gene modules (dictionary size in SMAF).
        lda1 (float): Regularization parameter for SMAF.
        lda2 (float): Regularization parameter for SMAF.
        dataset_dir (str): Directory where the dataset is stored.
    
    Returns:
        dict: Best measurement matrix and evaluation metrics.
    """
    try:
        # --- Main Simulation Pipeline ---
        
        # Load the AnnData object
        adata = sc.read_h5ad(adata_path)
        
        # Downsample and get filename
        downsampled_file = downsample_adata(adata, gene_set_size, num_cells, output_path=dataset_dir)
        
        # Load the downsampled data
        adata_subset = sc.read_h5ad(downsampled_file)
        
        # Use filename prefix for subsequent steps
        dataset_prefix = os.path.basename(downsampled_file).replace(".h5ad", "")
        
        # Split data and get list of generated files
        generated_files = split_and_save_data(adata_subset, dataset_prefix, dataset_dir)
        
        # Load split data
        train_data, validate_data, test_data = load_saved_data(dataset_prefix, dataset_dir)
        
        # Run SMAF
        U, W = smaf(train_data, num_modules, lda1=lda1, lda2=lda2, maxItr=20,
                    use_chol=False, donorm=True, mode=1, mink=0., doprint=True)
        
        # Remove zero-contribution modules
        nz = (U.sum(axis=0) > 0)
        U = U[:, nz]
        
        print("U dimensions =", U.shape)
        
        # Generate and evaluate measurement matrices
        best_coh_scores, Phi_coh = find_best_coherence_matrices(
            m=num_measurements, g=train_data.shape[0], 
            min_pools_per_gene=min_pools_per_gene, 
            max_pools_per_gene=max_pools_per_gene, 
            U=U,
            num_matrices=1, num_best=1,
        )
        
        best_rec_scores, Phi_best = find_best_reconstruction_matrices(
            Phi_coh, validate_data, U, sparsity, num_best=1
        )
        
        best_idx = np.argmax(best_rec_scores)
        best_matrix = Phi_best[best_idx]
        best_score = best_rec_scores[best_idx]

        return {
            "best_matrix": best_matrix,
            "best_score": best_score,
            "best_coherence_scores": best_coh_scores,
            "best_reconstruction_scores": best_rec_scores
        }

    finally:
        # --- Cleanup Section ---
        print("\nCleaning up temporary files...")

        def safe_delete(filepath):
            """Safely delete a file if it exists."""
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
                    print(f"Deleted: {filepath}")
            except Exception as e:
                print(f"Error deleting {filepath}: {e}")

        # Remove downsampled .h5ad file
        safe_delete(downsampled_file)

        # Remove generated .npy files
        for file_path in generated_files:
            safe_delete(file_path)

        print("Cleanup complete.\n")