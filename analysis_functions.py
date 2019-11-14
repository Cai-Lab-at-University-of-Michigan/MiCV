import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
sc.settings.autoshow = False

from app import cache

from helper_functions import *

def preprocess_data(session_ID, adata,
                    min_cells=2, min_genes=200, max_genes=10000,
                    target_sum=1e6, flavor="cell_ranger", 
                    n_top_genes=2000):

    adata = generate_adata_from_10X(session_ID)

    # do preprocessing
    # TODO: parametrize this in the webapp
    print("[STATUS] performing QC and normalizing data")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata


    # filter down to top 2000 highly-variable genes, using cell_ranger method
    # TODO: parametrize this in the webapp
    print("[STATUS] selecting highly variable genes")
    sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top_genes)
    sc.pl.highly_variable_genes(adata)
    adata = adata[:, adata.var['highly_variable']]

    # add some extra informational columns
    adata.obs["cell_ID"] = adata.obs.index
    adata.obs["cell_numeric_index"] = [i for i in range(0,len(adata.obs.index))]

    cache_adata(session_ID, adata)
    return adata

def do_PCA(session_ID, adata, n_comps=50, random_state=0):
    print("[STATUS] doing PCA")
    sc.tl.pca(adata, svd_solver="arpack", 
              n_comps=n_comps, random_state=random_state)
    new_adata = adata.copy()

    cache_adata(session_ID, new_adata)
    return new_adata

def do_neighborhood_graph(session_ID, adata, n_neighbors=20, random_state=0):
    print("[STATUS] finding neighbors")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    new_adata = adata.copy()

    cache_adata(session_ID, new_adata)
    return new_adata

def do_UMAP(session_ID, adata, random_state=0):
    print("[STATUS] doing a UMAP projection")
    new_adata = sc.tl.umap(adata, random_state=random_state, 
                           init_pos="spectral", copy=True)
    
    cache_adata(session_ID, new_adata)
    return new_adata

def do_clustering(session_ID, adata, resolution=0.5,
                  random_state=0, copy=True):
    print("[STATUS] performing clustering")
    sc.tl.leiden(adata, resolution=resolution, 
                 random_state=random_state)
    new_adata = adata.copy()
    new_adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])

    cache_adata(session_ID, new_adata)
    return new_adata