import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns


def scatter_pl(adata_in,type_specific=False,save_as=None,colorby=None,std=False):#optional plotting function
    for type in adata_in.obs['cell_type'].unique():
        adata=adata_in[adata_in.obs['cell_type']==type].copy()
        # Step 1: Define chromosome order map
        chrom_order_map = {
            **{str(i): i for i in range(1, 23)},  # chr1 - chr22
            "X": 23,
            "Y": 24,
            "M": 25,
            "MT": 25
        }

        # Step 2: Prepare var
        var = adata.var.copy()
        var['chrom'] = var['chromosome'].str.replace('chr', '', regex=False)
        var['chrom_order'] = var['chrom'].map(chrom_order_map).fillna(99).astype(int)
        var = var.dropna(subset=['start', 'end'])
        var['midpoint'] = ((var['start'] + var['end']) // 2).astype(int)

        # Step 3: Sort features by chromosome and midpoint
        var_sorted = var.sort_values(by=['chrom_order', 'midpoint'])
        #  Compute chromosome offsets
        chrom_sizes = var_sorted.groupby('chrom_order')['midpoint'].max().sort_index()
        chrom_offsets = chrom_sizes.cumsum().shift(fill_value=0)

        #  Map each gene to its genome-wide coordinate
        var_sorted['genome_coord'] = var_sorted.apply(
            lambda row: row['midpoint'] + chrom_offsets[row['chrom_order']],
            axis=1
        )

        gene_order = var_sorted.index.values
    
        # Step 4: Reorder X matrix and prepare dot plot data
        if type_specific==True:
            if std==True:
                X = adata[adata.obs['cell_type']==type, gene_order].layers['ctype log std norm']
            else:
                X = adata[adata.obs['cell_type']==type, gene_order].layers['ctype log norm']
        else:
            if std==True:
                X = adata[adata.obs['cell_type']==type, gene_order].layers['global log std norm']
            else:
                X = adata[adata.obs['cell_type']==type, gene_order].layers['global log norm']
        if hasattr(X, 'toarray'):
            X = X.toarray()

        n_cells = adata.shape[0]
        expression_values = X.flatten()
        genome_positions = np.repeat(var_sorted['genome_coord'].values, n_cells)
        cell_indices = np.tile(np.arange(n_cells), len(var_sorted))

        # Step 5: Plot

        if colorby!=None:
            #mapping colors to simulated_cnvs
            #Generate color map based on unique 'simulated_cnvs' values
            # Sanitize: convert CNV strings to just the first CNV per cell
            cnv_labels = adata.obs[colorby].astype(str).str.split(',').str[0].str.strip()

            # Create color map
            unique_labels = cnv_labels.unique()
            color_map = {label: plt.cm.tab20(i % 20) for i, label in enumerate(unique_labels)}

            # Assign colors per cell
            colors_for_cells = cnv_labels.map(color_map).values  # shape: (n_cells,)

            # Repeat the colors for each gene in each cell
            # Create the repeated list of colors (for each gene in each cell)
            colors = np.repeat(colors_for_cells, len(var_sorted))

        #Flatten expression values and genome positions as you did before
        n_cells = adata.shape[0]
        expression_values = X.flatten()  # Flattened expression values
        genome_positions = np.repeat(var_sorted['genome_coord'].values, n_cells)  # Flattened genome positions


        #actual plot
        plt.figure(figsize=(18, 6))
        if colorby!=None:
            plt.scatter(genome_positions, expression_values, s=0.2, alpha=0.2, rasterized=True,c=colors)
        else:
            plt.scatter(genome_positions, expression_values, s=0.2, alpha=0.2, rasterized=True)
        plt.xlabel("Genomic Position")
        plt.ylabel("Expression (Log2(exp/avg))")
        plt.ylim((-5,5))
        if type_specific==True:
            if std==True:
                plt.title(type + " Type Specific std Normalized Expression by Genomic Position")
            else:
                plt.title(type + " Type Specific Normalized Expression by Genomic Position")
        else:
            if std==True:
                plt.title(type + " Global std Normalized Expression by Genomic Position")
            else:
                plt.title(type + " Global Normalized Expression by Genomic Position")

        # Step 6: Chromosome separators and labels
        xticks = []
        xtick_labels = []

        for chrom_order in sorted(var_sorted['chrom_order'].unique()):
            chrom_df = var_sorted[var_sorted['chrom_order'] == chrom_order]
            mid = chrom_df['genome_coord'].median()
            right = chrom_df['genome_coord'].max()
            
            chrom_label = chrom_df['chrom'].iloc[0]
            xticks.append(mid)
            xtick_labels.append(f"chr{chrom_label}")

            plt.axvline(x=right, color='gray', linestyle='--', linewidth=0.5)

        # Apply x-ticks
        plt.xticks(xticks, xtick_labels, rotation=90, fontsize=8)
        plt.tight_layout()
        if save_as!=None:
            plt.savefig(save_as+type+'.png',dpi=300,bbox_inches='tight')
            plt.close()
            
        else:
            plt.show()
    return
