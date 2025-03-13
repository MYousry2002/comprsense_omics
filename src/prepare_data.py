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
from scipy.stats import spearmanr, pearsonr
import time
import os

THREADS = 10

# Get the absolute path of the project's root directory
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def downsample_adata(adata, gene_set_size, cell_count=10000, output_path=None):
    """
    Downsamples an AnnData object based on a specified gene set size.
    
    Parameters:
    - adata: AnnData object to be downsampled
    - gene_set_size: int, number of genes to include (500, 1000, or 5000)
    - cell_count: int, number of cells to select (default: 10,000)
    - output_path: str, directory for saving the downsampled file
    
    Returns:
    - str: Path to the downsampled AnnData file
    """
    assert gene_set_size in [500, 1000, 5000], "Gene set size must be 500, 1000, or 5000."

    # Default output path if not provided
    if output_path is None:
        output_path = os.path.join(ROOT_DIR, "dataset/pmotorcortex_mouse")

    # Ensure the output directory exists
    os.makedirs(output_path, exist_ok=True)

    # Load genes from file
    gene_file = os.path.join(ROOT_DIR, f"dataset/genes_{gene_set_size}.csv")
    
    if not os.path.exists(gene_file):
        raise FileNotFoundError(f"Gene file not found: {gene_file}")

    genes = pd.read_csv(gene_file, header=None)[1].tolist()

    # Normalize gene names
    adata.var_names = adata.var_names.str.upper().str.strip()
    genes = [gene.upper().strip() for gene in genes]

    # Match genes with existing ones in the dataset
    selected_genes = [gene for gene in genes if gene in adata.var_names]
    print(f"Matched genes for {gene_set_size} set: {len(selected_genes)}")

    # Ensure unique gene names
    adata.var_names_make_unique()

    # Prevent oversampling if cell_count is greater than available cells
    cell_count = min(cell_count, adata.n_obs)
    selected_cells = np.random.choice(adata.obs_names, size=cell_count, replace=False)

    # Subset the AnnData object
    adata_subset = adata[selected_cells, selected_genes]

    # Generate a unique filename
    timestamp = int(time.time())
    output_file = os.path.join(output_path, f"downsampled_{gene_set_size}genes_{os.getpid()}_{timestamp}.h5ad")

    # Save the new dataset
    adata_subset.write(output_file)
    print(f"Downsampled dataset saved to {output_file}")

    return output_file  # Return the filename


def split_and_save_data(adata, dataset_prefix, output_dir=None):
    """
    Splits `adata.X` into train, validation, and test subsets and saves them in the specified directory.
    
    Parameters:
    - adata: AnnData object, containing the dataset.
    - dataset_prefix: str, manually specified prefix for saved dataset files.
    - output_dir: str, directory to save the subsets.
    """
    if output_dir is None:
        output_dir = os.path.join(ROOT_DIR, "dataset/pmotorcortex/pmotorcortex_mouse")

    os.makedirs(output_dir, exist_ok=True)

    # Ensure `adata.X` is in a compatible format
    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X, dtype=np.float64)

    # Split dataset (50% Train, 25% Validate, 25% Test)
    num_cells = X.shape[0]
    train_idx, temp_idx = train_test_split(np.arange(num_cells), test_size=0.50, random_state=23)
    validate_idx, test_idx = train_test_split(temp_idx, test_size=0.50, random_state=23)

    # Extract data using indices
    train_data = X[train_idx]
    validate_data = X[validate_idx]
    test_data = X[test_idx]

    # Convert to dense arrays before saving
    train_data = np.asarray(train_data.todense() if sp.issparse(train_data) else train_data, dtype=np.float64)
    validate_data = np.asarray(validate_data.todense() if sp.issparse(validate_data) else validate_data, dtype=np.float64)
    test_data = np.asarray(test_data.todense() if sp.issparse(test_data) else test_data, dtype=np.float64)

    # Convert train_data to Fortran-contiguous format for SPAMS
    train_data = np.asfortranarray(train_data)

    # Save files
    file_paths = []
    for subset, data in zip(["train", "validate", "test"], [train_data, validate_data, test_data]):
        file_path = os.path.join(output_dir, f"{dataset_prefix}_{subset}_data.npy")
        np.save(file_path, data)
        file_paths.append(file_path)

    # Print dataset shapes
    print(f"Train data shape: {train_data.shape}")
    print(f"Validation data shape: {validate_data.shape}")
    print(f"Test data shape: {test_data.shape}")
    print(f"Datasets saved in {output_dir} with prefix '{dataset_prefix}'")

    return file_paths  # Return list of generated files


def load_saved_data(dataset_name, input_dir=None):
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
    if input_dir is None:
        input_dir = os.path.join(ROOT_DIR, "dataset/pmotorcortex/pmotorcortex_mouse")

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