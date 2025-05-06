import scanpy as sc
import numpy as np
import pandas as pd


def demo_heatmap(adata, cell_type_ann = 'cell_type', cell_type = str, chromosome_ann = 'chromosome', chromosome = str, cnv_ann = str,layers='counts'):
    # format match for chromosome index, some are int, some are like chr[int]
    if chromosome == 'X':
        chromosome = '23'
    if chromosome == 'Y':
        chromosome = '24'
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
    
    # gene-scaled heatmap makes the CNA footprint more apparent
    sc.pl.heatmap(adChr, adChr.var_names, groupby=cnv_ann, layer=layers, log=True, vmin=-1,vmax=1)
