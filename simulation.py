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
        dataset_dir (str): Directory where the dataset is stored.
    
    Returns:
        dict: Best measurement matrix and evaluation metrics.
    """
    
    # Load the AnnData object
    adata = sc.read_h5ad(adata_path)
    
    # Downsample the AnnData object
    adata_subset = downsample_adata(adata, gene_set_size, num_cells, output_path=dataset_dir)
    
    # Split the dataset into train, validation, and test sets
    split_and_save_data(adata_subset, "subset_data", dataset_dir)
    
    # Load the split data
    train_data, validate_data, test_data = load_saved_data("subset_data", dataset_dir)
    
    # Generate gene modules matrix U using SMAF
    U, W = smaf(train_data, num_modules, lda1=8, lda2=0.2, maxItr=20,
                use_chol=False, donorm=True, mode=1, mink=0., doprint=True)
    
    # Remove zero-contribution modules
    nz = (U.sum(axis=0) > 0)
    U = U[:, nz]
    
    print("U dimentions =", U.shape)
    
    # Generate and evaluate measurement matrices based on coherence
    best_coh_scores, Phi_coh = find_best_coherence_matrices(m=num_measurements, g=train_data.shape[0], 
        min_pools_per_gene=min_pools_per_gene, max_pools_per_gene=max_pools_per_gene, U=U,
                                                            num_matrices=5000, num_best=500,
    )
    
    # Find best measurement matrices based on reconstruction quality
    best_rec_scores, Phi_best = find_best_reconstruction_matrices(
        Phi_coh, validate_data, U, sparsity, num_best=50
    )
    
    # Select the best matrix based on reconstruction quality
    best_idx = np.argmax(best_rec_scores)
    best_matrix = Phi_best[best_idx]
    best_score = best_rec_scores[best_idx]
    
    return {
        "best_matrix": best_matrix,
        "best_score": best_score,
        "best_coherence_scores": best_coh_scores,
        "best_reconstruction_scores": best_rec_scores
    }
