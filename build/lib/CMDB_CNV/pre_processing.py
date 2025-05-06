import scanpy as sc
import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt
import scipy as sci
from scipy.sparse import issparse
import seaborn as sns


def pre_processing(ad10f, min_genes_on_chr = 5, keep_standard_chr_only = True):
    #preprocessing data
    #removed log norm step because we will be doing our own log normalizing
    #no genes were filtered out because i wanted maximum genome representation later
    ad10f.var['mt'] = ad10f.var_names.str.startswith('MT-')
    ribo_prefix = ("RPS","RPL")
    ad10f.var['ribo'] = ad10f.var_names.str.startswith(ribo_prefix)
    sc.pp.calculate_qc_metrics(ad10f, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
    
    adata = ad10f.copy()
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
    
    adClean = adClean[adClean.obs['pct_counts_mt']<20,:].copy()
    sc.pp.filter_cells(adClean, min_genes=500)
    sc.pp.filter_cells(adClean, max_counts=50000)
    sc.pp.filter_genes(adClean, min_cells=3)
    print('After filtering low-expression genes, ' + str(adClean.n_vars) + ' genes left.')

    adNorm = adClean.copy()
    adNorm.layers['counts'] = adNorm.X.copy()
    sc.pp.normalize_total(adNorm , target_sum=1e4)
    adNorm.layers['data'] = adNorm.X.copy()
    sc.pp.log1p(adNorm)
    adNorm.layers['lognorm'] = adNorm.X.copy()
    adNorm.X = adNorm.layers['data']
    return(adNorm.copy())
