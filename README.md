# CMDB_CNV: Copy Number Variation Analysis in scRNA-seq

`CMDB_CNV` is a Python package designed to detect and visualize copy number variations (CNVs) from single-cell RNA sequencing (scRNA-seq) data. It provides end-to-end tools for preprocessing, normalization, CNV neighborhood detection, and data visualization including heatmaps and scatter plots. 

---

## 📦 Installation

Clone the repository and install dependencies:

```bash
pip install git+https://github.com/Pooky-coding/CMDB_CNV.git
```

---

## 🚀 General Usage Overview

Import the main interface:

```python
from CMDB_CNV import CNV

adata = sc.read_h5ad('your h5ad file')
adata = CNV.simulate()                
adata = CNV.pre_processing           
adata = CNV.mean_norm(adata)
adata = CNV.mean_norm_std(adata)
adata = CNV.neighborhood(adata)
CNV.heatmap(adata)  
CNV.scatter_pl(adata)    
```

---

## 🧬 Functionality

### 🔹`simulate(adata, chrom, start, end, cn_state=2, fraction=0.1,cell_type=none,cell_type_column='cell_type',layer='counts',target_obs_col='simulated_cnvs_new')` 

Simulate a copy number alteration (CNA) in a specific genomic region for a fraction of cells in an AnnData object.

**Parameters:**

- `chrom` (*str or int*):
   The chromosome you want to simulate a CNV.
- `start` :
   The start(bp) of the CNV.
- `end` :
   The end (bp) of the CNV.
- `cn_state` (*int*, default=`2`):
   Integer representing the desired copy number.
- `fraction` (default=`0.1`):
   Fraction of eligible cells to apply the CNV to.
- `cell_type` (default=`'None'`):
   Restrict simulation to a specific cell type.
- `cell_type_column` (default=`'cell_type'`):
   Column in .obs with cell type info.
- `layer` (default=`'counts'`):
   name of the layer with counts (or use .X if not found).
- `target_obs_col` (default=`'simulated_cnvs_new'`):
   Column name in .obs to record simulated CNVs.

**Returns:**

- `AnnData` object (`anndata.AnnData`):
    Modified AnnData object with simulated CNA



------



### 🔹 `mean_norm(adata)`

Apply global and cell-type-specific mean normalization. Results are stored in various `adata.layers`.

**Normalization Strategies:**

1. **Global Mean Normalization (`adata.layers['log norm']`)**
   - Log-transforms the expression matrix: `log1p(expression)`
   - Subtracts the mean expression of each gene across all cells.
2. **Cell-Type Mean Normalization (`adata.layers['ctype log norm']`)**
   - Groups cells by `cell_type`.
   - Computes and subtracts gene-wise mean expression within each cell type.
   - Preserves expression variability within cell types while accounting for baseline differences.

**Parameters:**

- `adata` (*anndata.AnnData*):
   The input AnnData object. It must contain:
  - `.X`: Gene expression matrix (cells × genes)
  - `.obs['cell_type']`: Cell type annotations

**Returns:**

- `AnnData` object (`anndata.AnnData`):
   Modified `adata` with added layers:
  - `.layers['log norm']`: Globally mean-centered log-transformed expression
  - `.layers['ctype log norm']`: Cell-type mean-centered expression



------

### 🔹 `mean_norm_std(adata)`

Apply global and cell-type-specific mean normalization with ± 3 standard deviation. Results are stored in various `adata.layers`.

**Normalization Strategies:**

1. **Global Mean Normalization (`adata.layers['log norm']`)**
   - Log-transforms the expression matrix: `log1p(expression)`
   - Subtracts the mean expression of each gene across all cells.
2. **Cell-Type Mean Normalization (`adata.layers['ctype log norm']`)**
   - Groups cells by `cell_type`.
   - Computes and subtracts gene-wise mean expression within each cell type.
   - Preserves expression variability within cell types while accounting for baseline differences.

**Parameters:**

- `adata` (*anndata.AnnData*):
   The input AnnData object. It must contain:
  - `.X`: Gene expression matrix (cells × genes)
  - `.obs['cell_type']`: Cell type annotations

**Returns:**

- `AnnData` object (`anndata.AnnData`):
   Modified `adata` with added layers:
  - `.layers['log norm']`: Globally mean-centered log-transformed expression
  - `.layers['ctype log norm']`: Cell-type mean-centered expression

------


### 🔹 `pre_processing(adata)`

Clean and format AnnData object for CNV analysis.

**Processing Steps:**

- Removes genes on mitochondrial (`MT`) and ribosomal chromosomes (typically 'chrM', 'chrR', etc.)
- Sorts genes within each chromosome by start position.
- Computes a linear gene order index across all chromosomes.
- Validates that each gene has defined `chromosome`, `start`, and `end` fields.

**Parameters:**

`adata` (*anndata.AnnData*):
 The input AnnData object. It must contain:

- `.X`: Gene expression matrix (cells × genes)
- `.obs['cell_type']`: Cell type annotations
- `.var` with required columns:
  - `chromosome`: Chromosome identifiers (e.g., 'chr1', 'chr2', ..., 'chrX')
  - `start`: Gene start positions
  - `end`: Gene end positions

**Returns:**

- `AnnData` object (`anndata.AnnData`):
   Preprocessed `adata` with:
  - Filtered and ordered `.var` based on genomic coordinates.
  - Added `.var['gene_order']` field indicating gene positions along chromosomes.
  - Cleaned gene set for accurate CNV visualization and computation.



------



### 🔹 `neighborhood(adata, target_gene, window)`

Retrieve genes in a genomic window around a target gene to explore CNV neighborhoods. This is useful for visualizing and analyzing local CNV effects across contiguous gene regions.

**Parameters:**

- `adata` (*anndata.AnnData*):
   The input AnnData object. It must have:
  - `.var` with a `gene_order` field assigned during preprocessing
  - Proper genomic annotations (`chromosome`, `start`, `end`) in `.var`
- `target_gene` (*str*):
   The name (index in `adata.var_names`) of the gene of interest. This gene must exist in the dataset.
- `window` (*int*, default=`10`):
   Number of genes to include upstream and downstream of the target gene, resulting in a neighborhood of `2 × window + 1` genes.

**Returns:**

- `pandas.DataFrame`:
   A DataFrame containing the metadata for the neighborhood genes, including their:
  - Gene names (as index)
  - Chromosome, start, end
  - `gene_order` value

**Notes:**

- The function automatically handles edge cases (e.g., if the target gene is near the start or end of a chromosome).
- Make sure `CNV.pre_processing` has been applied to the dataset before using this function.



------



### 🔹 `heatmap(adata, cell_type='all', layer='ctype log norm')`

Draws a side-by-side heatmap comparing multiple cell types to demonstrate differences in expression patterns across chromosomes. This is useful for visually identifying potential CNV regions specific to certain cell types.

**Parameters:**

- `adata` (*anndata.AnnData*):
   The input AnnData object. It must include:
  - `.obs['cell_type']`: Cell type annotations
  - `.var['chromosome']` and `.var['gene_order']`: Genomic information
  - A valid expression `layer` (e.g., `'ctype log norm'`)
- `layer` (*str*, default=`'ctype log norm'`):
   Expression layer to visualize. Must exist in `adata.layers`.
- `cell_types` (*List[str]*, optional):
   List of specific cell types to include. If `None`, uses all unique cell types in `adata.obs['cell_type']`.
- `figsize` (*tuple*, default=`(14, 6)`):
   Size of the entire figure (width, height). Each heatmap will be adjusted within this space.
- `cmap` (*str*, default=`'RdBu_r'`):
   Color map for the heatmaps (e.g., red = amplification, blue = deletion).

**Returns:**

- `None`
   A heatmap plot is shown directly using `matplotlib.pyplot`.

**Notes:**

- Make sure to run `CNV.pre_processing` and `CNV.mean_norm` beforehand.
- This function is ideal for quick visual comparison of potential CNVs between cell populations.



------



### 🔹 `scatter_pl(adata)`

Creates a scatter plot showing the expression of a specific gene across different cell types. This visualization helps detect potential CNV signals, such as consistent over- or under-expression of a gene in specific cell populations.

**Parameters:**

- `adata` (*anndata.AnnData*):
   The input AnnData object. Must include:
  - `.obs['cell_type']`: Cell type annotations
  - Expression data in the specified `layer`
- `gene` (*str*):
   The name of the gene to visualize. Must be present in `adata.var_names`.
- `layer` (*str*, default=`'zscore'`):
   The expression layer to use for plotting (e.g., `'zscore'`, `'ctype log norm'`). This should be a normalized expression matrix.
- `figsize` (*tuple*, default=`(8, 5)`):
   Dimensions of the plot (width, height).
- `jitter` (*float*, default=`0.25`):
   Amount of horizontal jitter to apply to points to avoid overlap.

**Returns:**

- `None`
   A scatter plot is displayed directly using `matplotlib.pyplot`.



**Notes:**

- This plot is particularly useful for visualizing expression distribution within and across cell types.

- For best results, apply `CNV.mean_std_norm` before using this function to produce the `'zscore'` layer.

  

---

## 📁 Module Structure

| Function           | Description                                  |
| ------------------ | -------------------------------------------- |
| `simulate()`       | Simulates scRNA-seq data                     |
| `mean_norm()`      | Mean normalization                           |
| `mean_std_norm()`  | Z-score normalization                        |
| `pre_processing()` | Prepares the AnnData object                  |
| `neighborhood()`   | Extracts gene neighborhoods                  |
| `heatmap()`        | Creates chromosome-based expression heatmaps |
| `scatter_pl()`     | Plots gene expression in neighborhoods       |

---

## 📌 Data Requirements

The input `AnnData` object (`adata`) should include:

- `.X`: Expression matrix
- `.obs['cell_type']`: Cell labels
- `.var` with `chromosome`, `start`, and `end` columns

---

## 🧪 Dependencies

- `scanpy`
- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`

---

## 👨‍💻 Author

Developed by [Tye Chicha, Mingxuan Chi, Yingshan(Sam) Bi, Jiaxin Li].  
Open to feedback, issues, and contributions!
