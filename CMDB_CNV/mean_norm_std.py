def mean_norm_std(adata):
    """
    Applies both global and cell-type-specific mean ±3*std normalization to the input AnnData object.
    
    Outputs:
    - 'col_mean_std_norm' and 'global log std norm' layers for global normalization
    - 'ctype_mean_std_norm' and 'ctype log std norm' layers for cell-type-specific normalization
    """
    # Apply global normalization
    adata = global_mean_std_norm(adata)
    
    # Apply cell-type-specific normalization
    adata = ctype_mean_std_norm(adata)
    
    return adata

def global_mean_std_norm(adata):#generate global mean +- 3sd normalized counts
    Y = adata.X
    if not isinstance(Y, np.ndarray) and issparse(Y):
        Y = Y.toarray()  # Convert sparse matrix to dense if needed

    # Calculate the mean and standard deviation for each column
    col_means = np.mean(Y, axis=0)  # Mean of each column
    col_stds = np.std(Y, axis=0)    # Standard deviation of each column

    # Avoid division by zero for standard deviations
    col_stds[col_stds == 0] = 1e-10  # Prevent division by zero if std is zero
    
    # Compute the normalization factors:
    col_mean_plus_3std = col_means + 3 * col_stds
    col_mean_minus_3std = col_means - 3 * col_stds

    # Initialize normalized matrix
    Y_normalized = np.zeros_like(Y)

    # Normalize values based on the condition
    for i in range(Y.shape[1]):  # Iterate over each column (gene)
        for j in range(Y.shape[0]):  # Iterate over each row (cell)
            if Y[j, i] < col_means[i]:
                # If value is less than the column mean, divide by (mean - 3 * std)
                Y_normalized[j, i] = Y[j, i] / col_mean_minus_3std[i]
            else:
                # If value is greater than or equal to the column mean, divide by (mean + 3 * std)
                Y_normalized[j, i] = Y[j, i] / col_mean_plus_3std[i]

    # Optionally apply a log transformation (if needed)
    Y_log2_norm = np.log2(Y_normalized)  # Use +1 to avoid log(0) issues
    Y_log2_norm[np.isneginf(Y_log2_norm)]=0
    # Store the results in layers
    adata.layers['col_mean_std_norm'] = Y_normalized
    adata.layers['global log std norm'] = Y_log2_norm
    
    return adata


def ctype_mean_std_norm(adata):##generate cell type specific mean +- 3sd normalized counts
    adata.layers['ctype_mean_std_norm'] = adata.X.copy()  # Store original data if needed
    adata.layers['ctype log std norm'] = adata.X.copy()       # Store original data if needed
    
    n_cells, n_genes = adata.shape
    result_log2 = np.zeros((n_cells, n_genes))
    result_mean = np.zeros((n_cells, n_genes))
    
    for cell_type in adata.obs['cell_type'].unique():
        idx = adata.obs['cell_type'] == cell_type
        Y = adata[idx].X.copy()  # Get the data for the current cell type
        if not isinstance(Y, np.ndarray) and issparse(Y):
            Y = Y.toarray()

        # Calculate column-wise means and standard deviations
        col_means = np.mean(Y, axis=0)  # Mean of each column (gene)
        col_stds = np.std(Y, axis=0)    # Standard deviation of each column (gene)

        # Avoid division by zero for standard deviations
        col_stds[col_stds == 0] = 1e-10  # Prevent division by zero if std is zero

        # Compute the normalization factors
        col_mean_plus_3std = col_means + 3 * col_stds
        col_mean_minus_3std = col_means - 3 * col_stds

        # Initialize normalized matrix
        Y_normalized = np.zeros_like(Y)

        # Normalize values based on the condition
        for i in range(Y.shape[1]):  # Iterate over each gene (column)
            for j in range(Y.shape[0]):  # Iterate over each cell (row)
                if Y[j, i] < col_means[i]:
                    # If value is less than the column mean, divide by (mean - 3 * std)
                    Y_normalized[j, i] = Y[j, i] / col_mean_minus_3std[i]
                else:
                    # If value is greater than or equal to the column mean, divide by (mean + 3 * std)
                    Y_normalized[j, i] = Y[j, i] / col_mean_plus_3std[i]

        # Optionally apply a log transformation (if needed)
        Y_log2_norm = np.log2(Y_normalized)  # Use +1 to avoid log(0) issues
        Y_log2_norm[np.isneginf(Y_log2_norm)] = 0  # Replace -inf values with 0

        # Store the results for this cell type
        result_log2[idx.to_numpy(), :] = Y_log2_norm
        result_mean[idx.to_numpy(), :] = Y_normalized

    # Store the final results in the adata object layers
    adata.layers['ctype log std norm'] = result_log2
    adata.layers['ctype_mean_std_norm'] = result_mean

    return adata
