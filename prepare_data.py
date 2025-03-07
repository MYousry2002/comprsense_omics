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


def downsample_adata(adata, gene_set_size, cell_count=10000, output_path="dataset/cb_adult_mouse"):
    """
    Downsamples an AnnData object based on a specified gene set size.
    
    Parameters:
    - adata: AnnData object to be downsampled
    - gene_set_size: int, number of genes to include (500, 1000, or 5000)
    - cell_count: int, number of cells to select (default: 10,000)
    - output_path: str, base path for saving the downsampled file
    
    Returns:
    - AnnData object with selected genes and cells
    """
    assert gene_set_size in [500, 1000, 5000], "Gene set size must be 500, 1000, or 5000."
    
    # Load genes from file
    gene_file = f"./dataset/genes_{gene_set_size}.csv"
    genes = pd.read_csv(gene_file, header=None)[1].tolist()
    
    # Normalize gene names
    adata.var_names = adata.var_names.str.upper().str.strip()
    genes = [gene.upper().strip() for gene in genes]
    
    # Match genes with existing ones in the dataset
    selected_genes = [gene for gene in genes if gene in adata.var_names]
    print(f"Matched genes for {gene_set_size} set: {len(selected_genes)}")
    
    # Ensure unique gene names
    adata.var_names_make_unique()
    
    # Randomly select cells
    selected_cells = np.random.choice(adata.obs_names, size=cell_count, replace=False)
    
    # Subset the AnnData object
    adata_subset = adata[selected_cells, selected_genes]
    
    # Save the new dataset
    output_file = f"{output_path}_{gene_set_size}genes.h5ad"
    adata_subset.write(output_file)
    print(f"Downsampled dataset saved to {output_file}")
    
    return adata_subset


def split_and_save_data(adata, dataset_prefix, output_dir="./dataset/"):
    """
    Splits `adata.X` into train, validation, and test subsets and saves them in the specified directory.
    
    Parameters:
    - adata: AnnData object, containing the dataset.
    - dataset_prefix: str, manually specified prefix for saved dataset files.
    - output_dir: str, directory to save the subsets.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Ensure `adata.X` is in a compatible format
    X = adata.X
    
    # Convert sparse matrices to compressed row format
    if sp.issparse(X):  
        X = X.tocsr()
    else:
        X = np.asarray(X, dtype=np.float64)
    
    # Get total number of samples (cells)
    num_cells = X.shape[0]
    
    # Define split sizes (50% Train, 25% Validate, 25% Test)
    train_idx, temp_idx = train_test_split(np.arange(num_cells), test_size=0.50, random_state=23)
    validate_idx, test_idx = train_test_split(temp_idx, test_size=0.50, random_state=23)
    
    # Extract data using indices
    train_data = X[train_idx]
    validate_data = X[validate_idx]
    test_data = X[test_idx]
    
    # Convert all subsets to dense arrays before saving
    train_data = np.asarray(train_data.todense() if sp.issparse(train_data) else train_data, dtype=np.float64)
    validate_data = np.asarray(validate_data.todense() if sp.issparse(validate_data) else validate_data, dtype=np.float64)
    test_data = np.asarray(test_data.todense() if sp.issparse(test_data) else test_data, dtype=np.float64)

    # Convert train_data to Fortran-contiguous format for SPAMS
    train_data = np.asfortranarray(train_data)
    
    # Save subsets with manually inputted prefix
    np.save(os.path.join(output_dir, f"{dataset_prefix}_train_data.npy"), train_data)
    np.save(os.path.join(output_dir, f"{dataset_prefix}_validate_data.npy"), validate_data)
    np.save(os.path.join(output_dir, f"{dataset_prefix}_test_data.npy"), test_data)
    
    # Print dataset shapes
    print(f"Train data shape: {train_data.shape}")
    print(f"Validation data shape: {validate_data.shape}")
    print(f"Test data shape: {test_data.shape}")
    print(f"Datasets saved in {output_dir} as {dataset_prefix}_train_data.npy, {dataset_prefix}_validate_data.npy, and {dataset_prefix}_test_data.npy")



def load_saved_data(dataset_name, input_dir="./dataset/"):
    """
    Loads the train, validation, and test datasets from the specified directory using the given dataset name.
    
    Parameters:
    - dataset_name: str, name of the dataset (prefix used when saving)
    - input_dir: str, directory containing the saved subsets
    
    Returns:
    - train_data: numpy array
    - validate_data: numpy array
    - test_data: numpy array
    """
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Directory '{input_dir}' does not exist.")

    # Construct file paths using the dataset name
    train_data_path = os.path.join(input_dir, f"{dataset_name}_train_data.npy")
    validate_data_path = os.path.join(input_dir, f"{dataset_name}_validate_data.npy")
    test_data_path = os.path.join(input_dir, f"{dataset_name}_test_data.npy")

    # Check if files exist before loading
    for file_path in [train_data_path, validate_data_path, test_data_path]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

    # Load data
    train_data = np.load(train_data_path, allow_pickle=False)
    validate_data = np.load(validate_data_path, allow_pickle=False)
    test_data = np.load(test_data_path, allow_pickle=False)

    # Ensure all loaded data is in dense format
    validate_data = np.asarray(validate_data, dtype=np.float64)
    test_data = np.asarray(test_data, dtype=np.float64)
    
    train_data, validate_data, test_data = train_data.T, validate_data.T, test_data.T

    # Print dataset shapes
    print(f"Loaded dataset: {dataset_name}")
    print(f"Train data shape: {train_data.shape}")
    print(f"Validation data shape: {validate_data.shape}")
    print(f"Test data shape: {test_data.shape}")
    
    return train_data, validate_data, test_data