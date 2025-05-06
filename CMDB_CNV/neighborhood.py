import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

def get_valid_gene_posi(var, gene_names):
    invalid_genes = list()
    valid_genes = list()
    for gene in gene_names:
        start = var['start'][gene]
        end = var['end'][gene]
        # print(start, end)
        if np.isnan(start) or np.isnan(end):
            invalid_genes.append(gene)
        else:
            valid_genes.append(gene)
    return (invalid_genes, valid_genes)

def insert_nan_cols_numpy(df, nan_inserts):
    """
    Inserts multiple NaN columns into a DataFrame at specified column positions using NumPy for speed.
    
    Parameters:
        df (pd.DataFrame): Original DataFrame.
        nan_inserts (list of tuples): Each tuple is (num_nan_cols, loc).
    
    Returns:
        pd.DataFrame: New DataFrame with NaN columns inserted.
    """
    # Sort by insertion location
    nan_inserts_sorted = sorted(nan_inserts, key=lambda x: x[1])
    
    # Extract original values
    values = df.to_numpy()
    colnames = list(df.columns)
    n_rows = df.shape[0]
    
    new_values_parts = []
    new_colnames = []
    last_col = 0

    for num_nan_cols, loc in nan_inserts_sorted:
#         print(loc)
        # Append data before insertion point
        new_values_parts.append(values[:, last_col:loc])
        new_colnames.extend(colnames[last_col:loc])

        # Add NaN columns
        nan_block = np.full((n_rows, num_nan_cols), np.nan)
        new_values_parts.append(nan_block)
        new_colnames.extend([f'nan_col_{loc}_{i+1}' for i in range(num_nan_cols)])

        last_col = loc

    # Append remaining columns
    new_values_parts.append(values[:, last_col:])
    new_colnames.extend(colnames[last_col:])

    # Concatenate final values and rebuild DataFrame
    final_values = np.hstack(new_values_parts)
    return pd.DataFrame(final_values, columns=new_colnames, index=df.index)


def shifted_mean_and_consistency_window(df, window_length_by_unit, smooth = False):
    arr = df.to_numpy()
    n_rows, n_cols = arr.shape
    half_window = window_length_by_unit // 2

    # Prepare two empty DataFrames to store the results
    mean_arr = np.full_like(arr, np.nan, dtype=float)
    consistency_arr = np.full_like(arr, np.nan, dtype=float)

    for i in range(n_rows):  # For each cell
        row = arr[i, :]
        shifted_values = []

        # Create shifts in both directions
        for s in range(-half_window, half_window + 1):
            shifted = np.roll(row, s)
            if s < 0:
                shifted[:abs(s)] = np.nan
            elif s > 0:
                shifted[-s:] = np.nan
            shifted_values.append(shifted)

        shifted_stack = np.stack(shifted_values, axis=0)  # shape: (window_size, n_bins)
        
        # Center row is always at index `half_window`
        center = shifted_stack[half_window, :]
        neighbors = np.delete(shifted_stack, half_window, axis=0)

        # Calculate mean for each bin (standard smoothing)
        if smooth:
            mean_arr[i, :] = np.nanmean(shifted_stack, axis=0)

        # Calculate consistency: 1 - mean(|center - neighbor|), ignoring NaNs
        abs_diffs = np.abs(neighbors - center)
        consistency_arr[i, :] = 1 - np.nanmean(abs_diffs, axis=0)

    # Convert back to DataFrame and return both results
    consistency_df = pd.DataFrame(consistency_arr, index=df.index, columns=df.columns)
    del consistency_arr
    if smooth:
        mean_df = pd.DataFrame(mean_arr, index=df.index, columns=df.columns)
        del mean_arr
        return mean_df, consistency_df
    else:
        return consistency_df



def concatenate_subsets_back(target_matrix, subset_matrix):
    # Ensure df_matrix_Mono_nc is a subset of df_matrix
    if not set(subset_matrix.columns).issubset(set(target_matrix.columns)):
        raise ValueError("subset_matrix must be a subset of df_matrix")
    
    # Fill in the columns of df_matrix_Mono_nc into the full_matrix
    target_matrix[subset_matrix.columns] = subset_matrix

    return target_matrix





'''
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
main fucntion
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
'''
def neighborhood(adata, fc_matrix_layer = 'ctype log norm', window_length_by_unit = 10, unit = 'kb', smooth = False):
    '''
    Calculates the neighborhood consistency (NC) of gene expression across the genome for each cell
    by smoothing expression values with their genomic neighbors along each chromosome.

    This function:
    - Takes a gene expression matrix from a specified layer of an AnnData object.
    - Inserts NaN spacer columns between genes according to their genomic distance.
    - Applies a shifting window to smooth each gene’s expression by averaging it with its neighbors.
      Outputs a new layer named '[]_smoothed' in the AnnData object.
    - Applies a shifting window to calculate local similarity score of CNA fold changes.
      Outputs a new layer named 'neighborhood_consistency' in the AnnData object.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with gene metadata in `adata.var` including 'start', 'end', and 'chromosome'.
    fc_matrix_layer : str, optional
        The name of the layer in `adata.layers` that contains the gene expression matrix to be smoothed.
        Default is 'ctype log norm'.
    window_length_by_unit : int, optional
        The number of units (e.g. kb) to determine how many neighboring columns to include on each side.
        Must be an odd integer to ensure a center point exists.
        Default is 10.
    unit : str, optional
        Unit of genomic distance used to calculate spacing between genes.
        Options include: 'hb' (hectobase), 'kb' (kilobase), 'dkb' (decakilobase), 'hkb' (hectokilobase), 'mb' (megabase).
        Default is 'kb'.
    smooth : Boolean, optional
        whether to perform gene expressin smoothing or not

    Returns:
    --------
    adata : AnnData
        The input AnnData object with a new layer 'neighborhood_consistency' containing local CNA consistency and '[]_smoothed' containing the smoothed expression matrix.

    Notes:
    ------
    - Assumes `adata.var` has valid genomic coordinates.
    - Genes with invalid or missing genomic information are excluded.
    - Internally uses optimized numpy-based operations for memory efficiency.
    '''

    unit_dict = {'hb':100, 'kb':1000, 'dkb':10000, 'hkb':100000, 'mb':1000000}
    spacer_length = int(unit_dict[unit])

    cell_names = adata.obs_names
    gene_names = adata.var_names

    fc_matrix = adata.layers[fc_matrix_layer].copy()
    df_matrix = pd.DataFrame(fc_matrix, index=cell_names, columns=gene_names)
    # Create an empty DataFrame with the same shape as df_matrix, initialized with NaNs
    df_matrix_nc = pd.DataFrame(np.nan, index=df_matrix.index, columns=df_matrix.columns)
    if smooth:
        df_matrix_smooth = pd.DataFrame(np.nan, index=df_matrix.index, columns=df_matrix.columns)

    var = adata.var.copy()

    (invalid_genes, valid_genes) = get_valid_gene_posi(var, gene_names)
    var['midpoint'] = ((var['start'] + var['end']) // 2).astype(int)

    # calculate neighborhood consistency chr by chr
    chromosomes = var['chromosome'].unique().tolist()
    for chr_ in tqdm(chromosomes, desc='Calculating NC per chromosome'):
        print('Start processing ' + chr_ + ':')
        varMono = var[var['chromosome'] == chr_].copy()
        gene_names = varMono.index
        # get a sorted gene order according to start posi
        varMono_sorted = varMono.sort_values(by='start', ascending=True)
        gene_names_sorted = varMono_sorted.index
        # subset the matrix in to active chr as well as sort the genes in order of position
        df_matrix_Mono_sorted = df_matrix[gene_names_sorted]

        nan_inserts = list()
        for gene_index in range(len(gene_names_sorted)-1):
            gene = gene_names_sorted[gene_index]
            gene_next = gene_names_sorted[gene_index + 1]
            start = varMono_sorted.loc[gene, 'start']
            start_next = varMono_sorted.loc[gene_next, 'start']
            # distance between two neighbor genes
            distance = start_next - start
            num_nan_cols = int(distance // spacer_length)
            nan_inserts.append((num_nan_cols, gene_index))
        # print(nan_inserts)

        # insert cols of nans as spacer of 1kb between genes
        df_matrix_Mono_sorted_spacer = insert_nan_cols_numpy(df_matrix_Mono_sorted, nan_inserts)
        del df_matrix_Mono_sorted
        print('spacer done!')

        # duplicate and shift rows by windows and calculate means as neighborhoos consistency
        if smooth:
            (df_matrix_Mono_sorted_spacer_smooth, df_matrix_Mono_sorted_spacer_nc) = shifted_mean_and_consistency_window(df=df_matrix_Mono_sorted_spacer, window_length_by_unit=window_length_by_unit, smooth=smooth)
            print('Smooth done!')
            print('NC calculation done!')
        else:
            df_matrix_Mono_sorted_spacer_nc= shifted_mean_and_consistency_window(df=df_matrix_Mono_sorted_spacer, window_length_by_unit=window_length_by_unit, smooth=smooth)
            print('NC calculation done!')
        del df_matrix_Mono_sorted_spacer

        # save the neighborhood consistency df back to fc_matrix
        df_matrix_Mono_nc = df_matrix_Mono_sorted_spacer_nc[gene_names]
        del df_matrix_Mono_sorted_spacer_nc
        df_matrix_nc = concatenate_subsets_back(target_matrix=df_matrix_nc, subset_matrix=df_matrix_Mono_nc)
        if smooth:
            df_matrix_Mono_smooth = df_matrix_Mono_sorted_spacer_smooth[gene_names]
            del df_matrix_Mono_sorted_spacer_smooth
            df_matrix_smooth = concatenate_subsets_back(target_matrix=df_matrix_smooth, subset_matrix=df_matrix_Mono_smooth)
        print(chr_ + ' done!')
        print()

    layer_matrix_nc = df_matrix_nc.to_numpy()
    adata.layers[fc_matrix_layer + '_neighborhood_consistency'] = layer_matrix_nc
    if smooth:
        layer_matrix_smooth = df_matrix_smooth.to_numpy()
        adata.layers[fc_matrix_layer + '_smoothed'] = layer_matrix_smooth

    return adata
