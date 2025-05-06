# tests/test_CNV.py

from CMDB_CNV import CNV
import numpy as np
import anndata

def create_dummy_adata():
    X = np.random.rand(10, 5)
    obs = {"cell_type": ["A"] * 5 + ["B"] * 5}
    var = {"gene_symbols": [f"gene_{i}" for i in range(5)]}
    return anndata.AnnData(X=X, obs=obs, var=var)

def test_mean_std_norm():
    adata = create_dummy_adata()
    result = CNV.mean_std_norm(adata)
    assert 'col_mean_std_norm' in result.layers
    assert 'global log std norm' in result.layers
    assert 'ctype_mean_std_norm' in result.layers
    assert 'ctype log std norm' in result.layers

def test_simulate_runs():
    adata = create_dummy_adata()
    chrom = 'chr1'  # Example chromosome
    start = 1000     # Example start position
    end = 2000       # Example end position
    try:
        result = CNV.simulate(adata, chrom, start, end)
    except Exception:
        assert False, "simulate() raised an exception"

