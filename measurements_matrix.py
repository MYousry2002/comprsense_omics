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

from utilities import get_observations, compare_results



def random_double_balanced(m, g, max_pools_per_gene, min_pools_per_gene):
    """
    Generates a random measurement matrix (`phi`) with "double balanced" characteristics.
    
    This function ensures:
    - Each gene is assigned to `max_pools_per_gene` pools initially.
    - Genes with fewer than `min_pools_per_gene` assignments are reassigned.
    - The resulting matrix is normalized so that each row sums to 1.

    Parameters:
    - m : int
        Number of pools (rows in the matrix).
    - g : int
        Number of genes (columns in the matrix).
    - max_pools_per_gene : int
        Maximum number of pools a gene can be assigned to.
    - min_pools_per_gene : int
        Minimum number of pools a gene must be assigned to.

    Returns:
    - phi : np.ndarray
        A (m × g) measurement matrix where each gene is assigned to a set of pools.
    """

    # Initialize an empty measurement matrix (pools × genes)
    phi = np.zeros((m, g))

    # Randomly assign each gene to pools up to `max_pools_per_gene` times
    for i in range(max_pools_per_gene):
        idx = np.random.choice(g, g, replace=False)  # Shuffle gene indices
        idx = idx % m  # Ensure indices fit within the pool size (modulo operation)
        phi[idx, np.arange(g)] = 1  # Assign genes to random pools

    # Ensure each gene is assigned to at least `min_pools_per_gene` pools
    for i in np.where(phi.sum(0) < min_pools_per_gene)[0]:  # Identify under-assigned genes
        p = phi.sum(1).max() - phi.sum(1)  # Compute imbalance for each pool
        p[np.where(phi[:, i])[0]] = 0  # Exclude pools that already contain the gene

        # If there are pools available for reassignment
        if p.sum() > 0:
            p = p / p.sum()  # Normalize probability distribution
            num_to_assign = min((p > 0).sum(), int(min_pools_per_gene - phi[:, i].sum()))
            idx = np.random.choice(m, num_to_assign, replace=False, p=p)  # Select new pools
            phi[idx, i] = 1  # Assign gene to additional pools

    # Normalize the matrix so that each row sums to 1 (avoid bias in pooling)
    phi = (phi.T / phi.sum(1)).T

    return phi



def find_best_coherence_matrices(m, g, max_pools_per_gene, min_pools_per_gene, U, num_matrices=5000, num_best=500):
    """
    Generates a set of random measurement matrices and selects the `num_best` matrices
    with the lowest (best) 90th percentile coherence scores.

    Parameters:
        num_matrices (int): Total number of random matrices to generate.
        num_best (int): Number of best matrices to retain.
        num_pools (int): Number of pools (rows in the measurement matrix).
        num_features (int): Number of features (columns in the measurement matrix).
        num_min (int): Minimum number of pools assigned per feature.
        num_max (int): Maximum number of pools assigned per feature.
        U (np.ndarray): Feature transformation matrix for coherence computation.

    Returns:
        tuple: (best coherence scores, corresponding best measurement matrices)
    """
    best = np.ones(num_best)  # Initialize best coherence scores (worst possible score = 1)
    Phi_coh = [None for _ in best]  # Corresponding best measurement matrices

    for x in range(num_matrices):  
        if np.mod(x, num_matrices // 10) == 0:  # Print progress every 10% of iterations
            print(f"Iteration {x} of {num_matrices}")

        # Generate a new random measurement matrix
        phi = random_double_balanced(m, g, max_pools_per_gene, min_pools_per_gene)  

        # Compute the 90th percentile of the cosine distance between projected feature vectors
        coh_90 = np.percentile(1 - distance.pdist(phi.dot(U).T, 'cosine'), 90)

        # If the new `phi` has a better coherence score, replace the worst-performing matrix
        if coh_90 < best.max():  
            i = np.argmax(best)  # Find the index of the worst current coherence score
            best[i] = coh_90  # Replace the worst score with the new, better coherence score
            Phi_coh[i] = phi  # Store the corresponding measurement matrix
    
    return best, Phi_coh



def find_best_reconstruction_matrices(Phi_coh, validate_data, U, sparsity=0.02, num_best=50):
    """
    Selects the best measurement matrices based on their ability to recover original gene expression patterns.

    Parameters:
        Phi_coh (list): List of best measurement matrices from coherence selection.
        validate_data (np.ndarray): Original validation dataset.
        U (np.ndarray): Feature transformation matrix.
        sparsity (float): Sparsity constraint for sparse decoding (default: 0.02).
        num_best (int): Number of best matrices to retain.

    Returns:
        tuple: (best reconstruction scores, corresponding best measurement matrices)
    """
    best = np.zeros(num_best)  # Initialize best reconstruction scores (worst possible score = 0)
    Phi = [None for _ in best]  # Corresponding best measurement matrices

    for phi in Phi_coh:
        # Generate simulated observations by applying `phi` to the validation data
        y = get_observations(validate_data, phi, snr=5)  # Simulate pooled noisy measurements
        
        # Use sparse decoding (LASSO) to recover gene module activations
        w = sparse_decode(y, phi.dot(U), sparsity, method='lasso')  
        
        # Reconstruct gene expression using the estimated module activations
        x2 = U.dot(w)  # Approximate gene expression matrix from recovered module weights
        
        # Compare reconstructed expression (`x2`) with the original validation data (`validate_data`)
        r = compare_results(validate_data, x2)  
        
        # If the new measurement matrix `phi` produces a better reconstruction, update `best` and `Phi`
        if r[2] > best.min():  # If the new reconstruction score is better than the worst stored score
            i = np.argmin(best)  # Find the index of the worst-performing matrix
            best[i] = r[2]  # Replace the worst reconstruction score with the new, better score
            Phi[i] = phi  # Store the corresponding measurement matrix
    
    return best, Phi