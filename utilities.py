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



def random_phi_subsets_g(m, g, n, d_thresh=0.4):
    """
    Generates a random measurement matrix (Phi) with m pools and g genes while ensuring minimal correlation.
    
    Parameters:
        m (int): Number of pools (samples).
        g (int): Number of genes.
        n (tuple): Range (min, max) for the number of pools each gene is assigned to.
        d_thresh (float): Maximum allowable correlation threshold between gene assignments.
    
    Returns:
        np.ndarray: A (m x g) matrix where each column represents gene assignment to pools.
    """
    Phi = np.zeros((m, g))  # Initialize measurement matrix
    Phi[np.random.choice(m, np.random.randint(n[0], n[1]+1), replace=False), 0] = 1  # Assign first gene
    
    for i in range(1, g):  # Assign pools for remaining genes
        dmax = 1
        while dmax > d_thresh:
            p = np.zeros(m)
            p[np.random.choice(m, np.random.randint(n[0], n[1]+1), replace=False)] = 1  # Randomly assign pools
            dmax = 1 - distance.cdist(Phi[:, :i].T, [p], 'correlation').min()  # Ensure minimal correlation
        
        Phi[:, i] = p  # Assign gene to pools
    
    Phi = (Phi.T / Phi.sum(1)).T  # Normalize by row sum
    return Phi

def get_observations(X0, Phi, snr=5, return_noise=False):
    """
    Simulates observations by applying the measurement matrix (Phi) to the true signal (X0) with added noise.
    
    Parameters:
        X0 (np.ndarray): Original gene expression matrix (genes x samples).
        Phi (np.ndarray): Measurement matrix mapping genes to pools.
        snr (float): Signal-to-noise ratio.
        return_noise (bool): If True, return noise separately.
    
    Returns:
        np.ndarray: Noisy observed measurements (pools x samples).
        (Optional) np.ndarray: Noise matrix.
    """
    noise = np.array([np.random.randn(X0.shape[1]) for _ in range(X0.shape[0])])
    noise *= np.linalg.norm(X0) / np.linalg.norm(noise) / snr  # Scale noise level
    
    if return_noise:
        return Phi.dot(X0 + noise), noise
    else:
        return Phi.dot(X0 + noise)

def compare_distances(A, B, random_samples=[], s=200, pvalues=False):
    """
    Compares the Euclidean distances between corresponding columns of matrices A and B.
    
    Parameters:
        A (np.ndarray): First data matrix.
        B (np.ndarray): Second data matrix.
        random_samples (list, optional): Boolean mask of selected samples.
        s (int): Number of samples to randomly select if not provided.
        pvalues (bool): Whether to return p-values for the correlation tests.
    
    Returns:
        tuple: Pearson and Spearman correlation coefficients (and p-values if pvalues=True).
    """
    if len(random_samples) == 0:
        random_samples = np.zeros(A.shape[1], dtype=bool)
        random_samples[:min(s, A.shape[1])] = True  # Select `s` random samples
        np.random.shuffle(random_samples)
    
    dist_x = distance.pdist(A[:, random_samples].T, 'euclidean')
    dist_y = distance.pdist(B[:, random_samples].T, 'euclidean')
    
    pear = pearsonr(dist_x, dist_y)  # Pearson correlation
    spear = spearmanr(dist_x, dist_y)  # Spearman correlation
    
    if pvalues:
        return pear, spear  # Return correlations with p-values
    else:
        return pear[0], spear[0]  # Return only correlation coefficients

def correlations(A, B):
    """
    Computes correlation metrics between two matrices A and B.
    
    Parameters:
        A (np.ndarray): First matrix.
        B (np.ndarray): Second matrix.
    
    Returns:
        tuple: Overall correlation, Spearman correlation, per-gene correlation, per-sample correlation.
    """
    p = (1 - distance.correlation(A.flatten(), B.flatten()))  # Overall correlation
    spear = spearmanr(A.flatten(), B.flatten())  # Spearman correlation
    
    dist_genes = np.zeros(A.shape[0])  # Per-gene correlation
    for i in range(A.shape[0]):
        dist_genes[i] = 1 - distance.correlation(A[i], B[i])
    
    pg = np.average(dist_genes[np.isfinite(dist_genes)])  # Mean per-gene correlation
    
    dist_sample = np.zeros(A.shape[1])  # Per-sample correlation
    for i in range(A.shape[1]):
        dist_sample[i] = 1 - distance.correlation(A[:, i], B[:, i])
    
    ps = np.average(dist_sample[np.isfinite(dist_sample)])  # Mean per-sample correlation
    
    return p, spear[0], pg, ps

def compare_results(A, B):
    """
    Aggregates multiple comparison metrics between matrices A and B.
    
    Parameters:
        A (np.ndarray): First matrix.
        B (np.ndarray): Second matrix.
    
    Returns:
        list: Combined correlation and distance comparison results.
    """
    results = list(correlations(A, B))[:-1]  # Get correlation metrics (excluding `ps`)
    results += list(compare_distances(A, B))  # Compare Euclidean distances (gene-wise)
    results += list(compare_distances(A.T, B.T))  # Compare Euclidean distances (sample-wise)
    return results