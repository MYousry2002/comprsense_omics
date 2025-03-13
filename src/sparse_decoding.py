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


def sparse_decode(Y, D, k, worstFit=1., mink=0, method='omp', nonneg=False):
    """
    Sparse decoding method - obtain module activations from composite measurements.
    
    Parameters:
    - Y : Observed measurement data (pools × samples)
    - D : Dictionary (module patterns) used for reconstruction
    - k : Sparsity constraint (number of nonzero coefficients)
    - worstFit : Minimum required reconstruction accuracy (default = 1)
    - mink : Minimum sparsity allowed (default = 0)
    - method : 'omp' (Orthogonal Matching Pursuit) or 'lasso' (L1 regularization)
    - nonneg : If True, forces non-negative solutions (for LASSO only)
    
    Returns:
    - W : Estimated module activations (modules × samples)
    """
    
    if method == 'omp':
        while k > mink:
            W = spams.omp(np.asfortranarray(Y), np.asfortranarray(D), L=k, numThreads=THREADS)
            W = np.asarray(W.todense())

            fit = 1 - np.linalg.norm(Y - D.dot(W))**2 / np.linalg.norm(Y)**2  # Compute fit accuracy

            if fit < worstFit:
                break  # Stop if fit is too low
            else:
                k -= 1  # Decrease sparsity constraint
                
            
    elif method == 'lasso':
        Ynorm = np.linalg.norm(Y)**2 / Y.shape[1]  # Normalize Y
        W = spams.lasso(np.asfortranarray(Y), np.asfortranarray(D),
                        lambda1=k * Ynorm, mode=1, numThreads=THREADS, pos=nonneg)
        W = np.asarray(W.todense())
    
    return W

