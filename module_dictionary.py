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

from sparse_decoding import sparse_decode




def smaf(X, d, lda1, lda2, maxItr=10, UW=None, posW=False, posU=True, use_chol=False, 
         module_lower=1, activity_lower=1, donorm=False, mode=1, mink=5, U0=[], 
         U0_delta=0.1, doprint=False):
    """
    Sparse Module Activity Factorization (SMAF) to learn gene modules.
    
    Parameters:
    - X: np.array, shape (genes, cells), gene expression matrix. Note: dims opposite of anndata 
    - d: int, number of gene modules (dictionary size)
    - lda1: float, regularization term for activity sparsity
    - lda2: float, regularization term for module sparsity
    - maxItr: int, number of iterations for optimization
    - UW: tuple, (U, W) initial matrices, optional
    - posW: bool, enforce non-negativity on W
    - posU: bool, enforce non-negativity on U
    - use_chol: bool, use Cholesky decomposition for faster computation
    - module_lower: float, lower bound for module entropy
    - activity_lower: float, lower bound for activity entropy
    - donorm: bool, normalize U matrix
    - mode: int, optimization mode
    - mink: int, minimum sparsity constraint
    - U0: list, optional initialization for U
    - U0_delta: float, delta for projected gradient descent
    - doprint: bool, print progress
    
    Returns:
    - U, W: np.arrays, learned dictionary and activity matrices
    """
    # Initialize U and W matrices if not provided
    if UW is None:
        U, W = spams.nmf(np.asfortranarray(X), return_lasso=True, K=d, numThreads=THREADS)
        W = np.asarray(W.todense())  # Convert sparse W to dense format
    else:
        U, W = UW  # Use provided matrices
    
    # Compute initial reconstruction of X
    Xhat = U.dot(W)
    # Compute normalization factor for regularization
    Xnorm = np.linalg.norm(X) ** 2 / X.shape[1]
    
    # Iterate to optimize U and W
    for itr in range(maxItr):
        if mode == 1:
            # Solve for U using Lasso regression with sparsity regularization
            U = spams.lasso(np.asfortranarray(X.T), D=np.asfortranarray(W.T),
                            lambda1=lda2 * Xnorm, mode=1, numThreads=THREADS,
                            cholesky=use_chol, pos=posU)
            U = np.asarray(U.todense()).T  # Convert to dense and transpose
        elif mode == 2:
            # Optionally use projected gradient descent if U0 is provided
            if len(U0) > 0:
                U = projected_grad_desc(W.T, X.T, U.T, U0.T, lda2, U0_delta, maxItr=400)
                U = U.T  # Transpose back
            else:
                # Solve for U using Lasso regression with different lambda settings
                U = spams.lasso(np.asfortranarray(X.T), D=np.asfortranarray(W.T),
                                lambda1=lda2, lambda2=0.0, mode=2, numThreads=THREADS,
                                cholesky=use_chol, pos=posU)
                U = np.asarray(U.todense()).T  # Convert to dense and transpose
        
        # Normalize U if required
        if donorm:
            U = U / np.linalg.norm(U, axis=0)
            U[np.isnan(U)] = 0  # Replace NaN values with zero
        
        # Solve for W
        if mode == 1:
            wf = (1 - lda2)  # Worst-fit tolerance for sparsity
            W = sparse_decode(X, U, lda1, worstFit=wf, mink=mink)
        elif mode == 2:
            if len(U0) > 0:
                W = projected_grad_desc(U, X, W, [], lda1, 0., nonneg=posW, maxItr=400)
            else:
                W = spams.lasso(np.asfortranarray(X), D=np.asfortranarray(U),
                                lambda1=lda1, lambda2=1.0, mode=2, numThreads=THREADS,
                                cholesky=use_chol, pos=posW)
                W = np.asarray(W.todense())  # Convert to dense
        
        # Compute updated reconstruction of X
        Xhat = U.dot(W)
        
        # Compute module and activity sizes based on entropy
        module_size = np.average([np.exp(entropy(abs(u))) for u in U.T if u.sum() > 0])
        activity_size = np.average([np.exp(entropy(abs(w))) for w in W.T])
        
        # Print progress if required
        if doprint:
            print(distance.correlation(X.flatten(), Xhat.flatten()), module_size, activity_size, lda1, lda2)
        
        # Adjust sparsity parameters dynamically
        if module_size < module_lower:
            lda2 /= 2.  # Decrease sparsity regularization for U
        if activity_size < activity_lower:
            lda2 /= 2.  # Decrease sparsity regularization for W
    
    return U, W  # Return learned matrices