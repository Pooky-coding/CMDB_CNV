import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns


def simulate(
    adata,
    chrom,
    start,
    end,
    cn_state=2,
    fraction=0.1,
    cell_type=None,
    cell_type_column="cell_type",
    layer="counts",
    target_obs_col="simulated_cnvs_new"
):
    """
    Simulate a copy number alteration (CNA) in a specific genomic region for a fraction
    of cells in an AnnData object.

    Parameters:
    - adata: AnnData object
    - chrom: chromosome (string or int, e.g., 'X' or 1)
    - start: start position of the region (bp)
    - end: end position of the region (bp)
    - cn_state: integer representing the desired copy number (e.g., 0, 1, 3)
    - fraction: fraction of eligible cells to apply the CNV to
    - cell_type: restrict simulation to a specific cell type (str)
    - cell_type_column: column in .obs with cell type info
    - layer: name of the layer with counts (or use .X if not found)
    - target_obs_col: column name in .obs to record simulated CNVs

    Returns:
    - Modified AnnData object with simulated CNA
    """

    adata_new = adata.copy()
    adata_new = adata_new[:, ~adata_new.var['chromosome'].isna() & ~adata_new.var['start'].isna() & ~adata_new.var['end'].isna()]
    var = adata_new.var.dropna(subset=['chromosome', 'start', 'end']).copy()

    # Ensure coordinates are integers
    var['start'] = var['start'].astype(int)
    var['end'] = var['end'].astype(int)

    # Handle chromosome type matching
    if var['chromosome'].dtype.kind in {'i', 'u'}:
        chrom = int(chrom)
    else:
        chrom = str(chrom)

    # Print gene span on chromosome
    var_chr = var[var['chromosome'] == chrom]
    if var_chr.empty:
        print(f"❗ No genes found on chromosome {chrom}.")
        return adata_new

    start_min = var_chr['start'].min()
    end_max = var_chr['end'].max()
    print(f"Chromosome {chrom} contains genes from {start_min:,} to {end_max:,} bp")

    # Find genes in the region
    gene_mask = (
        (var['chromosome'] == chrom) &
        (var['start'] >= start) &
        (var['end'] <= end)
    )
    target_genes = var.index[gene_mask]
    print(target_genes)

    if len(target_genes) == 0:
        print("❗ No genes found in the specified region.")
        return adata_new

    # Subset the adata object to keep valid genes only
    adata_new._inplace_subset_var(var.index)

    # Get eligible cells
    if cell_type:
        eligible_cells = adata_new.obs.index[adata_new.obs[cell_type_column] == cell_type]
        if len(eligible_cells) == 0:
            print(f"⚠️ No cells found with type '{cell_type}'. Using all cells instead.")
            eligible_cells = adata_new.obs.index
    else:
        eligible_cells = adata_new.obs.index

    n_target = int(len(eligible_cells) * fraction)
    if n_target == 0:
        print("❗ No eligible cells selected — check fraction or cell type.")
        return adata_new

    # Randomly select cells to modify
    np.random.seed(42)
    selected_cells = np.random.choice(eligible_cells, n_target, replace=False)

    # Prepare to apply CNA
    multiplier = cn_state / 2.0  # Assume diploid baseline = 2
    data_matrix = adata_new.layers[layer] if layer in adata_new.layers else adata_new.X

    # Convert gene and cell names to integer indices
    gene_indices = [adata_new.var_names.get_loc(g) for g in target_genes]
    cell_indices = [adata_new.obs_names.get_loc(c) for c in selected_cells]

    if not hasattr(data_matrix, "multiply"):  # Dense
        data_matrix[np.ix_(cell_indices, gene_indices)] *= multiplier
    else:  # Sparse
        submatrix = data_matrix[np.ix_(cell_indices, gene_indices)].multiply(multiplier)
        data_matrix[np.ix_(cell_indices, gene_indices)] = submatrix

    # Record the simulated CNV in obs
    cna_label = f"{chrom}:{start}-{end} (CN {cn_state})"
    if target_obs_col not in adata_new.obs.columns:
        adata_new.obs[target_obs_col] = ""

    for cell in selected_cells:
        current = adata_new.obs.at[cell, target_obs_col]
        adata_new.obs.at[cell, target_obs_col] = (
            cna_label if current == "" else f"{current}, {cna_label}"
        )

    print(f"✅ Simulated CNA {cna_label} in {n_target} '{cell_type or 'ALL'}' cells across {len(target_genes)} genes.")
    return adata_new
