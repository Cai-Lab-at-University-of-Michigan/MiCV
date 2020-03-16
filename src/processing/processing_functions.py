import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
sc.settings.autoshow = False
import scanpy.external as sce
import palantir

from helper_functions import *

def preprocess_data(session_ID, adata,
                    min_cells=2, min_genes=200, max_genes=10000,
                    target_sum=1e6, flavor="cell_ranger", 
                    n_top_genes=2000):

    print("[DEBUG] adata: " + str(adata))

    # do preprocessing
    print("[STATUS] performing QC and normalizing data")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


    # filter down to top 2000 highly-variable genes, using cell_ranger method
    print("[STATUS] selecting highly variable genes")
    sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top_genes)

    # add some extra informational columns
    adata.obs["cell_ID"] = adata.obs.index
    adata.obs["cell_numeric_index"] = [i for i in range(0,len(adata.obs.index))]
    print(["DEBUG: caching adata after preprocessing."])
    cache_adata(session_ID, adata)

    # save the gene list for fast lookup
    gene_list = adata.var.index.tolist()
    gene_list = [str(x) for x in gene_list]
    gene_list = list(sorted(gene_list, key=str.lower))
    cache_gene_list(session_ID, gene_list)
    return adata

def do_PCA(session_ID, adata, n_comps=50, random_state=0):
    print("[STATUS] doing PCA")
    sc.tl.pca(adata, svd_solver="arpack", 
              n_comps=n_comps, random_state=random_state)
    new_adata = adata.copy()

    cache_adata(session_ID, new_adata)
    return new_adata

def do_neighborhood_graph(session_ID, adata, method="standard",
                          n_neighbors=20, random_state=0):
    print("[STATUS] finding neighbors")
    if ((method == "standard") or (not ("batch" in adata.obs))):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    elif (method == "bbknn"):
        sc.external.pp.bbknn(adata, batch_key="batch", approx=True)

    new_adata = adata.copy()

    cache_adata(session_ID, new_adata)
    return new_adata

def do_UMAP(session_ID, adata, n_dim_proj=2, random_state=0):
    print("[STATUS] doing a UMAP projection")
    sc.tl.umap(adata, random_state=random_state, 
                           init_pos="spectral", n_components=2, 
                           copy=False)
    adata_3D = sc.tl.umap(adata, random_state=random_state, 
                           init_pos="spectral", n_components=3, 
                           copy=True)
    adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
    cache_adata(session_ID, adata)
    return adata

def do_clustering(session_ID, adata, resolution=0.5,
                  random_state=0, copy=True):
    print("[STATUS] performing clustering")
    sc.tl.leiden(adata, resolution=resolution, 
                 random_state=random_state)
    new_adata = adata.copy()
    new_adata.obs["leiden_n"] = pd.to_numeric(adata.obs["leiden"])

    cache_adata(session_ID, new_adata)
    return new_adata