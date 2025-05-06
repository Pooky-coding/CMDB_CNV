

import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns


def heatmap(adata,cell_type='all',layer='ctype log norm'):
    #generates heatmap rows=cells, col=genes, color == expression. 1 plot per cell type x chromosome comination
    #has option to limit to cell type. always does 1 plot per chromosome.
    # Step 1: Prepare gene metadata
    var = adata.var.copy()
    var = var.dropna(subset=['start', 'end'])  # Remove missing start/end values
    var['chrom'] = var['chromosome'].str.replace('chr', '', regex=False)

    # Define chromosome order
    chrom_order_map = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "M": 25, "MT": 25}
    var['chrom_order'] = var['chrom'].map(chrom_order_map).fillna(99).astype(int)
    var['midpoint'] = ((var['start'] + var['end']) // 2).astype(int)
    if cell_type=='all':
        # Step 2: Loop over cell types and chromosomes
        for ctype in adata.obs['cell_type'].unique():
            valid_chroms = [chrom for chrom in var['chrom'].unique() if chrom in chrom_order_map]
            i=0
            for chrom in sorted(valid_chroms, key=lambda x: chrom_order_map[x]):
                # Select genes for this chromosome
                var_chrom = var[var['chrom'] == chrom].sort_values(by='midpoint')
                gene_order = var_chrom.index.values

                # Skip if no genes
                if len(gene_order) == 0:
                    continue

                # Subset adata to this cell type and chromosome
                toplot = adata[adata.obs['cell_type'] == ctype, gene_order]

                # Get layer data
                data = toplot.layers[layer]
                if hasattr(data, 'toarray'):
                    data = data.toarray()

                # Build DataFrame
                df = pd.DataFrame(data, index=toplot.obs_names, columns=toplot.var_names)
                heatmap_data_clean = df.replace([np.inf, -np.inf], np.nan).fillna(0)

                # Skip if empty
                if heatmap_data_clean.shape[1] == 0 or heatmap_data_clean.shape[0] == 0:
                    continue
                if len(chrom)<6:
                    # Plot heatmap
                    sns.clustermap(
                        heatmap_data_clean,
                        row_cluster=True,
                        col_cluster=False,  # genes already ordered
                        cmap='bwr',
                        yticklabels=False,
                        xticklabels=False,
                        figsize=(12, 8),
                        vmin=-4,
                        vmax=4
                    )
                    plt.title(f"Cell Type: {ctype} | Chromosome: chr{chrom}")
                    plt.show()
                i=i+1
    else:
        valid_chroms = [chrom for chrom in var['chrom'].unique() if chrom in chrom_order_map]
        i=0
        for chrom in sorted(valid_chroms, key=lambda x: chrom_order_map[x]):
            # Select genes for this chromosome
            var_chrom = var[var['chrom'] == chrom].sort_values(by='midpoint')
            gene_order = var_chrom.index.values

            # Skip if no genes
            if len(gene_order) == 0:
                continue

            # Subset adata to this cell type and chromosome
            toplot = adata[adata.obs['cell_type'] == cell_type, gene_order]

            # Get layer data
            data = toplot.layers[layer]
            if hasattr(data, 'toarray'):
                data = data.toarray()

            # Build DataFrame
            df = pd.DataFrame(data, index=toplot.obs_names, columns=toplot.var_names)
            heatmap_data_clean = df.replace([np.inf, -np.inf], np.nan).fillna(0)

            # Skip if empty
            if heatmap_data_clean.shape[1] == 0 or heatmap_data_clean.shape[0] == 0:
                continue
            if len(chrom)<6:
                # Plot heatmap
                sns.clustermap(
                    heatmap_data_clean,
                    row_cluster=True,
                    col_cluster=False,  # genes already ordered
                    cmap='bwr',
                    yticklabels=False,
                    xticklabels=False,
                    figsize=(12, 8),
                    vmin=-4,
                    vmax=4
                )
                plt.title(f"Cell Type: {cell_type} | Chromosome: chr{chrom}")
                plt.show()
            i=i+1
    return
