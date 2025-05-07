# Final project report

🤓 This report documents all the code used to generate the relevant results and figures for each task. Most of the analyses and interpretations of the figures are included in the main paper.



## Overview

With the development of single-cell RNA sequencing (scRNA-seq) technologies, researchers can now explore genetic and epigenetic regulations at single-cell resolution. However, we can only obtain gene expression information from scRNA-seq data, which limits our understanding of genetic regulations at the DNA level. Copy Number Variants (CNVs), a major type of genetic variation, play essential roles in tumor progression. There is an increasing interest in inferring CNVs from gene expression in scRNA-seq. Several methods, such as InferCNV , have been developed. However, these methods are often restricted by stringent assumptions or dependencies on specific reference cells.

We developed **CMDB_CNV**, a Python-based package designed to infer CNVs from scRNA-seq data without requiring any reference cell annotation. The output from CMDB_CNV includes different visual plots, making it a versatile tool for analyzing CNVs from scRNA-seq datasets



## ⭐️ Task 1: Create a method that infers CNAs from scRNA-seq data

Here we delevelop a python package  `CMDB_CNV` , which is a Python package designed to detect and visualize copy number variations (CNVs) from single-cell RNA sequencing (scRNA-seq) data. It provides end-to-end tools for preprocessing, normalization, CNV neighborhood detection, and data visualization including heatmaps and scatter plots.

All the APIs for the package can be found in README file. More explannation on the mechanism can be seen in our paper. 

Below is the example code for using the package

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

from CMDB_CNV import CNV

#using the benchmark dataset as an example
ad = sc.read_h5ad('PBMC_simulated_cnas_041025.h5ad')
ad.X=pbmc.layers['counts'].copy() # this line is only for the benchmark dataset

ad_CNV = ad.copy()
ad_CNV = CNV.preprocessing(ad_CNV) # preprocessing the data for detecting CNV
ad_CNV = CNV.mean_norm(ad_CNV) # calculating the mean 
ad_CNV = CNV.neighborhood(ad_CNV,smooth=True,window_length_by_unit=2,unit='hkb') # need this step for heatmap generating 
CNV.scatter_pl(ad_CNV)   
CNV.heatmap(ad_CNV)  
```

---



## ⭐️ Task 2A:  Assessment

To evaluate the performance of our method, we predicted CNAs using benchmark datasets with both our tool and the published tool, inferCNV. We also examined how read depth affects the performance of our method in comparison to inferCNV.

### Part 1: Detecting CNVs from benchmark data 

#### 💻Code

```python
# using CMDB_CNV to analyzing benchmark data 
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

from CMDB_CNV import CNV

#using the benchmark dataset as an example
ad = sc.read_h5ad('PBMC_simulated_cnas_041025.h5ad')
ad.X=pbmc.layers['counts'].copy() # this line is only for the benchmark dataset

ad_CNV = ad.copy()
ad_CNV = CNV.preprocessing(ad_CNV) # preprocessing the data for detecting CNV
ad_CNV = CNV.mean_norm(ad_CNV) # calculating the mean 
ad_CNV = CNV.neighborhood(ad_CNV,smooth=True,window_length_by_unit=2,unit='hkb') # need this step for heatmap generating 
CNV.scatter_pl(ad_CNV)   
```

#### 📊Output


![pbmc_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/CMDB_CNV/pbmc_CD14%20monocyte_chr6_scatter.png)

![pbmc_CD14 monocyte_chr22_scatter](figure%20for%20report/Task2A/CMDB_CNV/pbmc_CD14%20monocyte_chr22_scatter.png)

![pbmc_CD14 monocyte_chr23_scatter](figure%20for%20report/Task2A/CMDB_CNV/pbmc_CD14%20monocyte_chr23_scatter.png)


### Part 2: Explor the impact of read depth

#### 1.  Generating Low-depth datasets 

Below is a code we used to down-sampling 

```python
import leidenalg
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scipy.sparse as sp
import anndata as ad
import pySingleCellNet as cn
import random
from collections import defaultdict
from scipy.stats import pearsonr

from CMDB_CNV import CNV

adBM = sc.read_h5ad('PBMC_simulated_cnas_041025.h5ad')
adBMNorm = CNV.preprocessing(adBM, min_genes_on_chr = 100)
raw_counts = adBMNorm.layers['counts']
avg_counts = raw_counts.sum(axis=1).mean()

# a downsampling function to alter the depth of datasets
def alter_depth(adata, counts = int, replace = False):
    adata_out = adata.copy()
    n_cells = adata.n_obs
    sc.pp.downsample_counts(adata_out, total_counts=counts*n_cells, replace=replace)  
    # Adjust target_counts as needed
    
    return adata_out
  
adBMNorm_500 = alter_depth(adBMNorm, counts=500)
adBMNorm_1000 = alter_depth(adBMNorm, counts=1000)
adBMNorm_2000 = alter_depth(adBMNorm, counts=2000)

# We therefore genrated a list of low-depth datasets.
```

#### 2. Analyzing low-depth datasets with CMDB_CNV

##### 💻Code 

```python
adBMNorm_500 = CNV.mean_norm(adBMNorm_500)
adBMNorm_500 = CNV.neighborhood(adBMNorm_500,smooth=True,window_length_by_unit=2,unit='hkb')
CNV.scatter_pl(adBMNorm_500) 
CNV.heatmap(adBMNorm_500)  

adBMNorm_1000 = CNV.mean_norm(adBMNorm_1000)
adBMNorm_1000 = CNV.neighborhood(adBMNorm_10000,smooth=True,window_length_by_unit=2,unit='hkb')
CNV.scatter_pl(adBMNorm_1000)
CNV.heatmap(adBMNorm_1000)  

adBMNorm_2000 = CNV.mean_norm(adBMNorm_2000)
adBMNorm_2000 = CNV.neighborhood(adBMNorm_2000,smooth=True,window_length_by_unit=2,unit='hkb')
CNV.scatter_pl(adBMNorm_2000)
CNV.heatmap(adBMNorm_2000)  
```

##### 📊Output

- 500 reads per cell: 

![c500_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/CMDB_CNV/c500_CD14%20monocyte_chr6_scatter.png)

![c500_CD14 monocyte_chr22_scatter](figure%20for%20report/Task2A/CMDB_CNV/c500_CD14%20monocyte_chr22_scatter.png)

![c500_CD14 monocyte_chr23_scatter](figure%20for%20report/Task2A/CMDB_CNV/c500_CD14%20monocyte_chr23_scatter.png)


- 1000 reads per cell:
  
![c1000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/CMDB_CNV/c1000_CD14%20monocyte_chr6_scatter.png)

![c1000_CD14 monocyte_chr22_scatter](figure%20for%20report/Task2A/CMDB_CNV/c1000_CD14%20monocyte_chr22_scatter.png)

![c1000_CD14 monocyte_chr23_scatter](figure%20for%20report/Task2A/CMDB_CNV/c1000_CD14%20monocyte_chr23_scatter.png)


- 2000 reads per cell: 

![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/CMDB_CNV/c2000_CD14%20monocyte_chr6_scatter.png)

![c2000_CD14 monocyte_chr22_scatter](figure%20for%20report/Task2A/CMDB_CNV/c2000_CD14%20monocyte_chr22_scatter.png)

![c2000_CD14 monocyte_chr23_scatter](figure%20for%20report/Task2A/CMDB_CNV/c2000_CD14%20monocyte_chr23_scatter.png)


#### 3. Analyzing low-depth datasets with inferCNVs

##### 💻Code

```python
cnv.tl.infercnv(
    adBMNorm_500,
    reference_key='cell_type',
    window_size=50,
)
cnv.pl.chromosome_heatmap(adBMNorm_500, groupby=["cell_type"])
cnv.pl.chromosome_heatmap(adBMNorm_500, groupby=['simulated_cnvs'])

cnv.tl.infercnv(
    adBMNorm_1000,
    reference_key='cell_type',
    window_size=50,
)
cnv.pl.chromosome_heatmap(adBMNorm_1000, groupby=["cell_type"])
cnv.pl.chromosome_heatmap(adBMNorm_1000, groupby=['simulated_cnvs'])

cnv.tl.infercnv(
    adBMNorm_2000,
    reference_key='cell_type',
    window_size=50,
)
cnv.pl.chromosome_heatmap(adBMNorm_2000, groupby=["cell_type"])
cnv.pl.chromosome_heatmap(adBMNorm_2000, groupby=['simulated_cnvs'])
```

##### 📊Output

- 500 reads per cell

![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/inferCNV/500_infer1.jpg)

![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/inferCNV/500_infer2.jpg)


- 1000 reads per cell

![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/inferCNV/1000_infer1.jpg)

![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/inferCNV/1000_infer2.jpg)


- 2000 reads per cell

  ![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/inferCNV/2000_infer1.jpg)

  ![c2000_CD14 monocyte_chr6_scatter](figure%20for%20report/Task2A/inferCNV/2000_infer2.jpg)


## ⭐️ Task 2B: Augment assessment with better gold standard data

### Part 1: Generating gold standard datasets

To robustly evaluate our method, we developed a Python-based simulation framework to introduce synthetic copy number variations (CNVs) into a benchmark single-cell RNA-seq dataset. This framework is now wraped as a API in `CMDB_CNV`  as `CNV.simulate` (see details in README.md)

Below is the code for generateing the gold standard data

```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import re

from CMDB_CNV import CNV

#using the benchmark dataset as an base dataset
ad = sc.read_h5ad('PBMC_simulated_cnas_041025.h5ad')

# Start with a clean copy of your base dataset
adata_base = ad.copy()

# Define CNA configurations
cna_configs = [
    # (chrom, start, end, cn_state, frequency, size_label, freq_label)
    ('7', 3_000_000, 6_000_000, 0, 0.1, 'small', 'low'), #3MB
    ('8', 10_000_000, 13_000_000, 1, 0.5, 'small', 'medium'),
    ('9', 20_000_000, 23_000_000, 4, 0.9, 'small', 'high'),

    ('X', 2_000_000, 12_000_000, 4, 0.1, 'medium', 'low'), # 10MB
    ('16', 6_000_000, 16_000_000, 0, 0.5, 'medium', 'medium'),
    ('22', 20_000_630, 30_000_630, 1, 0.9, 'medium', 'high'),

    ('1', 3_000_000, 33_000_000, 1, 0.1, 'large', 'low'), # 30MB
    ('13', 10_000_000, 43_000_000, 4, 0.5, 'large', 'medium'),
    ('20', 20_000_000, 50_000_000, 0, 0.9, 'large', 'high'),
]

# Define cell types to rotate through
cell_types = ['B cell', 'CD14 monocyte', 'CD4 T cell', 'CD8 T cell']

# Split CNVs into groups by frequency
low_cnas = [c for c in cna_configs if c[4] == 0.1]
medium_cnas = [c for c in cna_configs if c[4] == 0.5]
high_cnas = [c for c in cna_configs if c[4] == 0.9]

# Helper function to apply a list of CNAs to a copy of the base data
def simulate_and_save(cna_list, filename, adata_base, cell_types):
    adata_simulated = adata_base.copy()

    for i, (chrom, start, end, cn_state, freq, size_label, freq_label) in enumerate(cna_list):
        cell_type = cell_types[i % len(cell_types)]
        print(f"\n➡️ Simulating {size_label} CNA ({freq_label} freq): Chr{chrom}:{start}-{end}, CN={cn_state}, Cell type: {cell_type}")
        adata_simulated = CNV.simulate(
            adata=adata_simulated,
            chrom=chrom,
            start=start,
            end=end,
            cn_state=cn_state,
            fraction=freq,
            cell_type=cell_type,
            cell_type_column="cell_type",
            layer="counts",
            target_obs_col="simulated_cnvs_new"
        )

    # Save the simulated dataset
    adata_simulated.write(filename)
    print(f"💾 Saved: {filename}")

# Apply simulations and save each set
simulate_and_save(low_cnas, "simulated_cnas_low_freq.h5ad", adata_base, cell_types)
simulate_and_save(medium_cnas, "simulated_cnas_medium_freq.h5ad", adata_base, cell_types)
simulate_and_save(high_cnas, "simulated_cnas_high_freq.h5ad", adata_base, cell_types)

print("\n✅ All simulations complete and saved.")    
```

---

The code generated 3 different files. Each files contains 3 simulated CNVs with different size

![23431746555852_.pic](/Users/aurora/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/6fb8cf30fd4f4c5b3a1f5986804702d8/Message/MessageTemp/9e20f478899dc29eb19741386f9343c8/Image/23431746555852_.pic.jpg)

### Part 2: Detecting CNVs with CMDB_CNV

#### 💻 Code

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

from CMDB_CNV import CNV

#import the datasets
low_freq=sc.read_h5ad('simulated_cnas_low_freq.h5ad')
low_freq.X=low_freq.layers['counts'].copy()
med_freq=sc.read_h5ad('simulated_cnas_medium_freq.h5ad')
med_freq.X=med_freq.layers['counts'].copy()
high_freq=sc.read_h5ad('simulated_cnas_high_freq.h5ad')
high_freq.X=high_freq.layers['counts'].copy()

low_freq=CNV.preprocessing(low_freq)
med_freq=CNV.preprocessing(med_freq)
high_freq=CNV.preprocessing(high_freq)

low_freq=CNV.mean_norm(low_freq)
med_freq=CNV.mean_norm(med_freq)
high_freq=CNV.mean_norm(high_freq)

low_freq=CNV.neighborhood(low_freq)
med_freq=CNV.neighborhood(med_freq)
high_freq=CNV.neighborhood(high_freq)

CNV.scatter_pl(low_freq)
CNV.scatter_pl(med_freq)
CNV.scatter_pl(high_freq)

```

---

#### 📊Output

The output here are scatter plots for the chromosom with CNVs in the corelated cell types

##### 1: Detecting CNVs with low frequency and different size 

![lowlow_B cell_chr7_scatter](figure%20for%20report/Task2B/CMDB_CNV/low_B%20cell_chr7_scatter.png)

![low_CD4 T cell_chr1_scatter](figure%20for%20report/Task2B/CMDB_CNV/low_CD4%20T%20cell_chr1_scatter.png)

![low_CD14 monocyte_chr23_scatter](figure%20for%20report/Task2B/CMDB_CNV/low_CD14%20monocyte_chr23_scatter.png)


##### 2: Detecting CNVs with medium frequency and different size 

![med_B cell_chr8_scatter](figure%20for%20report/Task2B/CMDB_CNV/med_B%20cell_chr8_scatter.png)

![med_CD4 T cell_chr13_scatter](figure%20for%20report/Task2B/CMDB_CNV/med_CD4%20T%20cell_chr13_scatter.png)

![med_CD14 monocyte_chr16_scatter](figure%20for%20report/Task2B/CMDB_CNV/med_CD14%20monocyte_chr16_scatter.png)


##### 3: Detecting CNVs with high frequency and different size 

![high_B cell_chr9_scatter](figure%20for%20report/Task2B/CMDB_CNV/high_B%20cell_chr9_scatter.png)

![high_CD4 T cell_chr20_scatter](figure%20for%20report/Task2B/CMDB_CNV/high_CD4%20T%20cell_chr20_scatter.png)

![high_CD14 monocyte_chr22_scatter](figure%20for%20report/Task2B/CMDB_CNV/high_CD14%20monocyte_chr22_scatter.png)


### 3. Using inferCNV to detect CNVs

To compare the performance, we also used inferCNV to find CNVs in our gold standard datasets

#### 💻 Code

```python
# import required packages
import leidenalg
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scipy.sparse as sp
import anndata as ad
import pySingleCellNet as cn
import random
from collections import defaultdict
from scipy.stats import pearsonr
import infercnvpy as cnv
import CMDB_CNV as CNV

plt.rcParams['figure.dpi'] = 300

# Set the global seed
seed = 1135809
# Python built-in random module
random.seed(seed)
# NumPy
np.random.seed(seed)
```

```python
#here is all the function we used whenever we use inferCNV
def cluster_plot(adNorm, n_pcs = None, n_neighbors = 20):
    adCluster = adNorm.copy()

    # kNN cell-cell distances
    sc.pp.neighbors(adCluster, n_neighbors = n_neighbors, n_pcs = n_pcs)

    # leiden clustering
    # for finer clustering, after multiple tests, i chose 0.12 for resolution. in this way, the cells are clustered into fine enough clusters for detecting sub types, while it is not too fine that many clusters representing one same cell type.
    sc.tl.leiden(adCluster, resolution = 0.12)
    # analyze intercluster similarities, for better agreement between embedding and clustering
    sc.tl.paga(adCluster)
    sc.pl.paga(adCluster, plot = False)
    # run umap with paga optimization
    sc.tl.umap(adCluster, 0.25, init_pos = 'paga', n_components=3)
    # visualize clustering
    sc.pl.umap(adCluster, color = ['leiden'], alpha = 0.75, s = 15, legend_loc = 'on data')
    return adCluster


def visualize(ad, cluster_key, color_key, signature_gene_dict, title = 'UMAP'):
    # visualize
    ad.obs[color_key] = ad.obs[color_key].astype('category')
    sc.tl.dendrogram(ad, color_key)
    fig, (ax_1, ax_2) = plt.subplots(2, 1, figsize = (10, 14), gridspec_kw = {'hspace':0.5})
    ax_1_dict = sc.pl.umap(ad, color = [color_key], alpha = 1, s = 10, legend_loc = 'on data', ax = ax_1, show = False, legend_fontsize = 8, add_outline=False, outline_color =('white', 'black'), title = title)
    # set label white with black outline for clearer visualization
    for text in ax_1.texts:
        text.set_color('white')
        text.set_fontweight('bold')  # Bold the font
        text.set_path_effects([
            path_effects.Stroke(linewidth=1, foreground='black'),
            path_effects.Normal()  # Normal text
        ])
    ax_2_dict = sc.pl.dotplot(ad, signature_gene_dict, color_key, dendrogram = True,ax = ax_2, show = False)
    plt.show()

def sub_clustering(ad, cluster_list, reso = 0.1):
    leiden_key = 'leiden'
    for cluster in cluster_list:
        sc.tl.leiden(ad, resolution = reso, restrict_to = [leiden_key, [cluster]], key_added = leiden_key + "_R" + cluster)
        leiden_key = leiden_key + "_R" + cluster
    return (ad, leiden_key)


def annotating(ad, cell_dict, list_remove, column_key_ref, column_key_target = 'cell_type_ann'):
    # remove doublets and debris
    ad_ann = ad.copy()
    # ~ to negate the mask of isin()
    ad_ann = ad_ann[~ad_ann.obs[column_key_ref].isin(list_remove)].copy()

    # annotate the cell types
    for cell in ad_ann.obs.index:
        cluster =  ad_ann.obs.loc[cell, column_key_ref] # for each cell get its cluster number
        ad_ann.obs.loc[cell, column_key_target] = cell_dict[cluster] # find the cell_type in the dict and save to column "cell_type"
    print('done!')
    
    return ad_ann

def run_monochr_celltype(adata, cell_type_ann = 'cell_type', cell_type = str, chromosome_ann = 'chromosome', chromosome = str, ref_key = str, ref_cat = None, ax = None, group_by = str, show = True):
    # format match for chromosome index, some are int, some are like chr[int]
    if chromosome == 'X':
        chromosome = '23'
    if chromosome == 'Y':
        chromosomr = '24'
    list_chromosome_index = ['chr'+str(chromosome), chromosome]
    try:
        list_chromosome_index.append(int(chromosome))
    except ValueError:
        pass

    try:
        list_chromosome_index.append(float(chromosome))
    except ValueError:
        pass    
    
    adMono = adata[adata.obs[cell_type_ann] == cell_type].copy()
    adChr = adMono[:,adMono.var[chromosome_ann].isin(list_chromosome_index)].copy()
    
    if ref_cat != None:
        cnv.tl.infercnv(
        adChr,
        reference_key=ref_key,
        reference_cat=ref_cat,
        window_size=10)
    else:
        cnv.tl.infercnv(
        adChr,
        reference_key=ref_key,
        window_size=10)
  
    if ax == None:
        cnv.pl.chromosome_heatmap(adChr, groupby=[group_by])
    else:
        cnv.pl.chromosome_heatmap(adChr, groupby=[group_by], ax = ax, show = show)
```

```python
# Below is the working code for using inferCNV to detect CNVs in gold standard data
def run_sample(ad, ref_key = 'cell_type', ref_cat = None):
    adNorm = scpipe(ad)
    ad_cnv = adNorm.copy()
    
    cnv.tl.infercnv(
    ad_cnv,
    reference_key=ref_key,
    reference_cat=ref_cat,
    window_size=50
    )
    
    cnv.pl.chromosome_heatmap(ad_cnv, groupby=['simulated_cnvs_new'])
    cnv.pl.chromosome_heatmap(ad_cnv, groupby=['cell_type'])
    
    cnv.tl.pca(ad_cnv, n_comps = 100)
    cnv.pp.neighbors(ad_cnv, n_neighbors = 5, n_pcs = 100)
    cnv.tl.leiden(ad_cnv)
    
    cnv.tl.umap(ad_cnv)
    cnv.tl.cnv_score(ad_cnv)
    cnv.pl.chromosome_heatmap(ad_cnv, groupby=['cnv_leiden'])
    
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(11, 11))
    ax4.axis("off")
    ax6.axis('off')
    cnv.pl.umap(
        ad_cnv,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(ad_cnv, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(ad_cnv, color="cell_type", ax=ax3, show=False)
    cnv.pl.umap(ad_cnv, color='simulated_cnvs_new', ax = ax5)
    
    return ad_cnv
  

adHigh = sc.read_h5ad('data/simulated_cnas_high_freq.h5ad')
adHigh_cnv = run_sample(adHigh)

adMed = sc.read_h5ad('data/simulated_cnas_medium_freq.h5ad')
adMed_cnv = run_sample(adMed)

adLow = sc.read_h5ad('data/simulated_cnas_low_freq.h5ad')
adLow_cnv = run_sample(adLow)
```

#### 📊Output

##### 1: Detecting CNVs with low frequency and different size

![low_1](figure%20for%20report/Task2B/inferCNV/low_1.jpg)

![low_2](figure%20for%20report/Task2B/inferCNV/low_2.jpg)

![low_3](figure%20for%20report/Task2B/inferCNV/low_3.jpg)

![low_4](figure%20for%20report/Task2B/inferCNV/low_4.jpg)


##### 2: Detecting CNVs with medium frequency and different size

![medium_1](figure%20for%20report/Task2B/inferCNV/medium_1.jpg)

![medium_2](figure%20for%20report/Task2B/inferCNV/medium_2.jpg)

![medium_3](figure%20for%20report/Task2B/inferCNV/medium_3.jpg)

![medium_4](figure%20for%20report/Task2B/inferCNV/medium_4.jpg)


##### 3: Detecting CNVs with high frequency and different size

![high_1](figure%20for%20report/Task2B/inferCNV/high_1.jpg)

![high_2](figure%20for%20report/Task2B/inferCNV/high_2.jpg)

![high_3](figure%20for%20report/Task2B/inferCNV/high_3.jpg)

![high_4](figure%20for%20report/Task2B/inferCNV/high_4.jpg)


## ⭐️ Task3: Measure CNVs in PSCs

To really analyze our package's performance, we applied our package to 3 different established scRNA-seq datasets.

Below is the code for this task pre-processing steps for each raw data 

### Part 0：Function for pre-processing the raw data

```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import re
from pybiomart import Server
warnings.filterwarnings('ignore')
plt.rcParams['figure.dpi'] = 300

def report_missing_chromosomal_info(adata, columns=['chromosome', 'start', 'end']):
    """
    Reports the number and percentage of genes with missing chromosomal information.

    Parameters:
    - adata: AnnData object with chromosomal info in .var
    - columns: List of column names to check for missing values (default: ['chromosome', 'start', 'end'])

    Returns:
    - Dictionary with counts and percentage
    """
    num_missing = adata.var[columns].isna().any(axis=1).sum()
    total_genes = adata.var.shape[0]
    percent_missing = num_missing / total_genes
    print(f"🧬 Genes with missing chromosomal info: {num_missing} out of {total_genes}")
    print(f"📉 That's {percent_missing:.2%} of all genes.")
    return adata
def add_chromosomal_info_biomart_ensemble_ID(adata, gene_column='ensembl_ids'):
    """
    Maps ensemble ids  to chromosomal information using Ensembl BioMart.

    Parameters:
    - adata: AnnData object
    - gene_column: Name of the column in adata.var that holds Ensembl gene IDs (default: 'ensemble_ids')

    Returns:
    - AnnData object with new columns in `.var`: 'chromosome_name', 'start_position', 'end_position'
    """
    # Connect to Ensembl BioMart
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    
    # Query BioMart for gene location information
    query = dataset.query(attributes=["ensembl_gene_id", "chromosome_name", "start_position", "end_position"])
    chromosomal_info = pd.DataFrame(query)
    
    # Rename columns for consistency
    chromosomal_info.columns = ["ensembl_gene_id", "chromosome_name", "start_position", "end_position"]
    
    # Drop duplicates in case of multiple mappings
    chromosomal_info = chromosomal_info.drop_duplicates(subset=["ensembl_gene_id"])

    # Ensure gene IDs are str
    adata.var['ensembl_ids'] = adata.var[gene_column].astype(str)
    chromosomal_info['ensembl_gene_id'] = chromosomal_info['ensembl_gene_id'].astype(str)

    # Map values using ensembl_gene_id as lookup
    adata.var['chromosome'] = adata.var['ensembl_ids'].map(chromosomal_info.set_index('ensembl_gene_id')['chromosome_name'])
    adata.var['start'] = adata.var['ensembl_ids'].map(chromosomal_info.set_index('ensembl_gene_id')['start_position'])
    adata.var['end'] = adata.var['ensembl_ids'].map(chromosomal_info.set_index('ensembl_gene_id')['end_position'])

    return adata
def norm_clean_data(adata, min_genes_on_chr = 5, keep_standard_chr_only = True):
    # Keep genes with non-null and non-empty chromosome, start, and end
    valid = (
        adata.var['chromosome'].notna() &
        adata.var['start'].notna() &
        adata.var['end'].notna() &
        (adata.var['chromosome'].astype(str).str.strip() != '')
    )
    adata = adata[:, valid].copy()
    
    # Make sure start and end are numeric
    adata.var['start'] = pd.to_numeric(adata.var['start'], errors='coerce')
    adata.var['end'] = pd.to_numeric(adata.var['end'], errors='coerce')
    adata = adata[:, adata.var['start'].notna() & adata.var['end'].notna()].copy()
    

    # remove all genes that are in a chromosome with less that 5 genes on it
    chr_counts = adata.var['chromosome'].value_counts()
    valid_chromosomes = chr_counts[chr_counts >= min_genes_on_chr].index
    
    # plus extract genes on standard chromosomes
    if keep_standard_chr_only:
        standard_chr = {str(i) for i in range(1, 23)} | {'X', 'Y'}
        valid_chromosomes = [chr_ for chr_ in valid_chromosomes if chr_ in standard_chr]

    adClean = adata[:, adata.var['chromosome'].isin(valid_chromosomes)].copy()
    print('After filtering invalid chromosomes, ' + str(adClean.n_vars) + ' genes left.')
    
    # make the genome annotation formated
    adClean.var["chromosome"] = adClean.var["chromosome"].astype(str).apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    adClean.var["chromosome"] = adClean.var["chromosome"].replace({"chrX": "chr23", "chrY": "chr24"})
    adClean.var["start"] = adClean.var["start"].astype(int)
    adClean.var["end"] = adClean.var["end"].astype(int)
    
    adNorm = adClean.copy()
    
    # remove genes expressed in too few cells, representing non-informative genes
    sc.pp.filter_genes(adNorm, min_cells = 3)
    print('After filtering non-informative genes, ' + str(adNorm.n_vars) + ' genes left.')
    
    # save raw counts in 'counts' layer
    adNorm.layers['counts'] = adNorm.X.copy()
    
    # normalize data + save in 'data' layer
    sc.pp.normalize_total(adNorm, target_sum = 1e4)
    adNorm.layers['data'] = adNorm.X.copy()

    # save scaled data in 'data_scaled' layer
    sc.pp.log1p(adNorm)
    adNorm.layers['data_scaled'] = adNorm.X.copy()
    
    print(adNorm.layers)
    sc.pp.highly_variable_genes(adNorm, min_mean=0.0125, max_mean=6, min_disp=0.25) # default flavor is seurat

    sc.tl.pca(adNorm, use_highly_variable=True)
    sc.pl.pca_variance_ratio(adNorm, 50)
    
    adNorm.X = adNorm.layers['data']

    print('Preprocessing done!')
    print('There are ' + str(adNorm.n_obs) + ' cells in this dataset.')
    return adNorm
```

### Part 1: scRNA-seq data processing

#### Dataset 1: 2D-gastruloids 

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

adata_gastruloids = sc.read_10x_h5("/Users/aurora/Desktop/computational SCB/final/1-2D-gastruloids.h5")
adata_g = adata_gastruloids.copy()
adata_g = add_chromosomal_info_biomart_ensemble_ID(adata_g,gene_column='gene_ids')
adata_g = report_missing_chromosomal_info(adata_g)
adata_g = norm_clean_data(adata_g)
n_neighbors = 20
n_pcs = 10
sc.pp.neighbors(adata_g, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.leiden(adata_g,.1)
sc.tl.umap(adata_g, 0.25)
sc.pl.umap(adata_g,color=['leiden'], alpha=.75, s=15, legend_loc='on data')
```

```python
cell_dict = {' Advanced Mesoderm':['0'],'Nascent Mesoderm': ['1'], 'Amniotic Ectoderm': ['2'],'Posterior Primitive Streak': ['3'],'Anterior Primitive Streak': ['4'],'Definitive Endoderm': ['5']}
marker_genes_dict = {'Posterior Primitive Streak':['HAND1', 'DOK4', 'AFAP1L2'],'Anterior Primitive Streak': ['TPM1', 'IGFBP7', 'S100A11'], 'Amniotic Ectoderm': ['PHC1', 'CD24', 'SFRP2'],'Advanced Mesoderm': ['APLNR', 'HAS2', 'ITGA5'],'Nascent Mesoderm': ['LMNB1', 'FAM169A', 'PDHB'],'Definitive Endoderm': ['SOX17', 'CDH2', 'CXCR4']}

new_obs_name = 'cell_type'
adata_g.obs[new_obs_name] = np.nan

for i in cell_dict.keys():
    ind = pd.Series(adata_g.obs.leiden).isin(cell_dict[i])
    adata_g.obs.loc[ind,new_obs_name] = i

adata_g.obs['cell_type'] = adata_g.obs['cell_type'].astype("category")
sc.tl.dendrogram(adata_g, "cell_type")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), gridspec_kw={'wspace':0.4})
ax1_dict = sc.pl.umap(adata_g,color=['cell_type'], alpha=.75, s=10, legend_loc='on data', ax=ax1, show=False, title='1_cell type composition')
ax2_dict = sc.pl.dotplot(adata_g, marker_genes_dict, 'cell_type', dendrogram=True,ax=ax2, show=False)
plt.show()
```

![task3-2D-gastruloids](figure%20for%20report/Task3/processing/task3-2D-gastruloids.jpg)

#### Dataset 2: ipsc_neurons

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

adata_neurons = sc.read_10x_h5("/Users/aurora/Desktop/computational SCB/final/iPSC_neurons.h5")
ad_neurons = adata_neurons.copy()
ad_neurons = add_chromosomal_info_biomart_ensemble_ID(ad_neurons,gene_column='gene_ids')
ad_neurons = report_missing_chromosomal_info(ad_neurons)
adNorm_neurons = norm_clean_data(ad_neurons)
n_neighbors = 20
n_pcs = 10
sc.pp.neighbors(adNorm_neurons, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.leiden(adNorm_neurons,.1)
sc.tl.umap(adNorm_neurons, 0.5)
sc.pl.umap(adNorm_neurons,color=['leiden'], alpha=.75, s=15, legend_loc='on data')
```

```python
cell_dict = {'neurons_0':['0'],'neurons_1': ['1'], 'neurons_2': ['2']}

new_obs_name = 'cell_type'
adNorm_neurons.obs[new_obs_name] = np.nan

for i in cell_dict.keys():
    ind = pd.Series(adNorm_neurons.obs.leiden).isin(cell_dict[i])
    adNorm_neurons.obs.loc[ind,new_obs_name] = i

adNorm_neurons.obs['cell_type'] = adNorm_neurons.obs['cell_type'].astype("category")
sc.tl.dendrogram(adNorm_neurons, "cell_type")
sc.pl.umap(adNorm_neurons,color=['cell_type'], alpha=.75, s=10, legend_loc='on data', show=False, title='2_cell type composition')
plt.show()
```

![task3-neurons](figure%20for%20report/Task3/processing/task3-neurons.jpg)

#### Dataset 3:  iPSC-derived cardiomyocytes

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

adata_HSC = sc.read("/Users/aurora/Desktop/computational SCB/final/3-iPSC -> HSC/iPSC_to_HSC.h5ad")
print(adata_HSC)
ad_HSC = adata_HSC.copy()
ad_HSC = add_chromosomal_info_biomart_ensemble_ID(ad_HSC,gene_column='ensembl_ids')
ad_HSC = report_missing_chromosomal_info(ad_HSC)
adNorm_HSC = norm_clean_data(ad_HSC)
n_neighbors = 20
n_pcs = 30
sc.pp.neighbors(adNorm_HSC, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.leiden(adNorm_HSC,.1)
sc.tl.umap(adNorm_HSC, 0.4)
sc.pl.umap(adNorm_HSC,color=['leiden'], alpha=.75, s=15, legend_loc='on data')
cell_dict = {'Cardiomyocyte progenitors':['0'],'Endothelial cells': ['1'], 'Cardiac fibroblasts': ['2']}
marker_genes_dict = {'Cardiomyocyte progenitors':['KRT8', 'MEST', 'KRT18'],'Endothelial cells': ['GNG11','ESAM','IGFBP4'], 
                     'Cardiac fibroblasts': ['MALAT1', 'CDH13', 'COL4A2']}

new_obs_name = 'cell_type'
adata_g.obs[new_obs_name] = np.nan

for i in cell_dict.keys():
    ind = pd.Series(adNorm_HSC.obs.leiden).isin(cell_dict[i])
    adNorm_HSC.obs.loc[ind,new_obs_name] = i

adNorm_HSC.obs['cell_type'] = adNorm_HSC.obs['cell_type'].astype("category")
sc.tl.dendrogram(adNorm_HSC, "cell_type")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), gridspec_kw={'wspace':0.4})
ax1_dict = sc.pl.umap(adNorm_HSC,color=['cell_type'], alpha=.75, s=10, legend_loc='on data', ax=ax1, show=False, title='3_cell type composition')
ax2_dict = sc.pl.dotplot(adNorm_HSC, marker_genes_dict, 'cell_type', dendrogram=True,ax=ax2, show=False)
plt.show()
```

![task3-cadio](figure%20for%20report/Task3/processing/task3-cadio.jpg)

#### Dataset 4：kidney organoids

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

adata_kidney = sc.read("/Users/aurora/Desktop/computational SCB/final/2-kidney/kidney.h5ad")
print(adata_kidney)
adata1 = adata_kidney.copy()
adata1 = add_chromosomal_info_biomart_ensemble_ID(adata1)
adata1 = report_missing_chromosomal_info(adata1)
adNorm_kidney = norm_clean_data(adata1)
n_neighbors = 20
n_pcs = 10
sc.pp.neighbors(adNorm_kidney, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.leiden(adNorm_kidney,.1)
sc.tl.umap(adNorm_kidney, 0.25)
sc.pl.umap(adNorm_kidney,color=['leiden'], alpha=.75, s=15, legend_loc='on data')
cell_dict = {'kidney organoids_0 ':['0'],'kidney organoids_1': ['1'], 'kidney organoids_2': ['2'],'kidney organoids_3': ['3'],'kidney organoids_4': ['4'],'kidney organoids_5': ['5'], 'kidney organoids_6': ['6']}

new_obs_name = 'cell_type'
adNorm_kidney.obs[new_obs_name] = np.nan

for i in cell_dict.keys():
    ind = pd.Series(adNorm_kidney.obs.leiden).isin(cell_dict[i])
    adNorm_kidney.obs.loc[ind,new_obs_name] = i

adNorm_kidney.obs['cell_type'] = adNorm_kidney.obs['cell_type'].astype("category")
sc.tl.dendrogram(adNorm_kidney, "cell_type")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), gridspec_kw={'wspace':0.4})
ax1_dict = sc.pl.umap(adNorm_kidney,color=['cell_type'], alpha=.75, s=10, legend_loc='on data', ax=ax1, show=False, title='2_cell type composition')
plt.show()
```

![task3-kidney](figure%20for%20report/Task3/processing/task3-kidney.jpg)

### Part 2:  Using CMDB_CNV to detect CNVs

#### 💻Code

```python
import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns

from CMDB_CNV import CNV

#import the datasets
gastruloids = sc.read_h5ad('2D-gastruloids.h5ad')
neurons = sc.read_h5ad('ips_neurons.h5ad')
cardio = sc.read_h5ad('iPSC-derived cardiomyocytes.h5ad')
kidney = sc.read_h5ad('kidney organoids.h5ad')

gastruloids = CNV.preprocessing(gastruloids)
neurons = CNV.preprocessing(neurons)
cardio = CNV.preprocessing(cardio)
kidney = CNV.preprocessing(kidney)

gastruloids = CNV.mean_norm(gastruloids)
neurons = CNV.mean_norm(neurons)
cardio = CNV.mean_norm(cardio)
kidney = CNV.mean_norm(kidney)

gastruloids = CNV.neighborhood(gastruloids)
neurons = CNV.neighborhood(neurons)
cardio = CNV.neighborhood(cardio)
kidney = CNV.neighborhood(kidney)


CNV.heatmap(gastruloids)
CNV.heatmap(neurons)
CNV.heatmap(cardio)
CNV.heatmap(kidney)
```

---

#### 📊Output 


##### 3.1 CMDB_CNV with ips_neurons
![neurons_S](figure%20for%20report/Task3/CMDB_CNV/neurons_s.jpg)
![neurons_H](figure%20for%20report/Task3/CMDB_CNV/neurons_h.jpg)

##### 3.2 CMDB_CNV with kidney organoids
![kidney_S](figure%20for%20report/Task3/CMDB_CNV/kidney_S.jpg)
![kidney_H](figure%20for%20report/Task3/CMDB_CNV/kidney_h.jpg)

### Part 3 Using inferCNV to detect CNVs for comparizion

Here we also used inferCNV to detect CNVs in the datasets we processed. We did that to analyze the performance of our package since we are not sure if any of the datasets really have CNVs

#### 3.1 InferCNV with ipsc_neurons

```python
import leidenalg
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scipy.sparse as sp
import anndata as ad
import pySingleCellNet as cn
import random
from collections import defaultdict
from scipy.stats import pearsonr
import infercnvpy as cnv
from CMDB_CNV import CNV

plt.rcParams['figure.dpi'] = 300

# Set the global seed
seed = 1135809
# Python built-in random module
random.seed(seed)
# NumPy
np.random.seed(seed)

adipsneuronNorm = sc.read_h5ad('data/ips_neurons.h5ad')
adipsneuron_cnv = adipsneuronNorm.copy()
cnv.tl.infercnv(
    adipsneuron_cnv,
    reference_key='cell_type',
    window_size=50,
)
cnv.pl.chromosome_heatmap(adipsneuron_cnv, groupby=["cell_type"])
cnv.pl.chromosome_heatmap_summary(adipsneuron_cnv, groupby="cell_type")
```

![task3-inferCNV_neurons_1](figure%20for%20report/Task3/inferCNV/task3-inferCNV_neurons_1.jpg)

![task3-inferCNV_neurons_2](figure%20for%20report/Task3/inferCNV/task3-inferCNV_neurons_2.jpg)



#### 3.2 InferCNV with kidney organoids

```python
import leidenalg
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scipy.sparse as sp
import anndata as ad
import pySingleCellNet as cn
import random
from collections import defaultdict
from scipy.stats import pearsonr
import infercnvpy as cnv
from CMDB_CNV import CNV

plt.rcParams['figure.dpi'] = 300

# Set the global seed
seed = 1135809
# Python built-in random module
random.seed(seed)
# NumPy
np.random.seed(seed)

adkidneyorgNorm = sc.read_h5ad('data/kidney organoids.h5ad')
adkidneyorg_cnv = adkidneyorgNorm.copy()
cnv.tl.infercnv(
    adkidneyorg_cnv,
    reference_key='cell_type',
    window_size=50,
)
cnv.pl.chromosome_heatmap(adkidneyorg_cnv, groupby=["cell_type"])
cnv.pl.chromosome_heatmap_summary(adkidneyorg_cnv, groupby="cell_type")
```

---

![task3-inferCNV_kidney](figure%20for%20report/Task3/inferCNV/task3-inferCNV_kidney.jpg)

![task3-inferCNV_kidney_2](figure%20for%20report/Task3/inferCNV/task3-inferCNV_kidney_2.jpg)


#### 3.3 InferCNV with  iPSC-derived cardiomyocytes

```python
import leidenalg
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scipy.sparse as sp
import anndata as ad
import pySingleCellNet as cn
import random
from collections import defaultdict
from scipy.stats import pearsonr
import infercnvpy as cnv
from CMDB_CNV import CNV

plt.rcParams['figure.dpi'] = 300

# Set the global seed
seed = 1135809
# Python built-in random module
random.seed(seed)
# NumPy
np.random.seed(seed)

adipscardiNorm = sc.read_h5ad('data/iPSC-derived cardiomyocytes.h5ad')
adipscardi_cnv = adipscardiNorm.copy()
cnv.tl.infercnv(
    adipscardi_cnv,
    reference_key='cell_type',
    window_size=50,
)
cnv.pl.chromosome_heatmap(adipscardi_cnv, groupby=["cell_type"])
cnv.pl.chromosome_heatmap_summary(adipscardi_cnv, groupby="cell_type")
```

---
![task3-inferCNV_cardi_1](figure%20for%20report/Task3/inferCNV/task3-inferCNV_cardi_1.jpg)

![task3-inferCNV_cardi_2](figure%20for%20report/Task3/inferCNV/task3-inferCNV_cardi_2.jpg)
