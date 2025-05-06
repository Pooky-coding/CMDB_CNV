import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns


def mean_norm(adata):
    """
    Applies both global and cell type-specific mean normalization
    to the input AnnData object and stores the results in layers.

    Returns:
        Updated AnnData object with additional normalized layers.
    """
    # Run global normalization
    adata = global_mean_norm(adata)
    
    # Run cell type-specific normalization
    adata = ctype_mean_norm(adata)
    
    return adata

def global_mean_norm(adata):#generates global mean normalized counts
    Y = adata.X.copy()
    if not isinstance(Y, np.ndarray):
        Y = Y.toarray()

    # Calculate column-wise non-zero means
    nonzero_counts = (Y != 0).sum(axis=0)     # number of non-zero values per column
    nonzero_sums = Y.sum(axis=0)              # sum of all values per column
    nonzero_counts[nonzero_counts == 0] = 1e-10  # prevent division by zero
    col_means_nonzero = nonzero_sums / nonzero_counts
    col_means_nonzero[col_means_nonzero==0]=1e-10 # safeguard against columns that have no nonzero values so their mean==0
    # Normalize each element by its column's non-zero mean
    Y_normalized = Y / col_means_nonzero
    Y_log2_norm=np.log2(Y_normalized)
    Y_log2_norm[np.isneginf(Y_log2_norm)]=-4
    # Replace or store
    adata.layers['col_mean_norm'] = Y_normalized#adata.X = X_normalized  # or 
    adata.layers['global log norm']=Y_log2_norm
    return(adata)



def ctype_mean_norm(adata):#generates cell type specific normalizzed counts
    adata.layers['ctype_col_mean_norm'] = adata.X.copy()#adata.X = X_normalized  # or 
    adata.layers['ctype log norm']=adata.X.copy()
    n_cells, n_genes = adata.shape
    result_log2 = np.zeros((n_cells, n_genes))
    result_mean = np.zeros((n_cells, n_genes))
    for type in adata.obs['cell_type'].unique():
        idx = adata.obs['cell_type'] == type
        Y = adata[idx].X
        if not isinstance(Y, np.ndarray):
            Y = Y.toarray()

        # Calculate column-wise non-zero means
        nonzero_counts = (Y != 0).sum(axis=0)     # number of non-zero values per column
        nonzero_sums = Y.sum(axis=0)              # sum of all values per column
        nonzero_counts[nonzero_counts == 0] = 1e-10  # prevent division by zero
        col_means_nonzero = nonzero_sums / nonzero_counts
        col_means_nonzero[col_means_nonzero==0]=1e-10 # safeguard against columns that have no nonzero values so their mean==0
        # Normalize each element by its column's non-zero mean
        Y_normalized = Y / col_means_nonzero
        Y_log2_norm=np.log2(Y_normalized)
        Y_log2_norm[np.isneginf(Y_log2_norm)]=-4
        # store
        result_log2[idx.to_numpy(), :] = Y_log2_norm
        result_mean[idx.to_numpy(), :] = Y_normalized
    adata.layers['ctype log norm'] = result_log2
    adata.layers['ctype_col_mean_norm'] = result_mean
    return(adata)


