import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
sc.settings.autoshow = False
import scanpy.external as sce
import palantir

#from app import cache

from helper_functions import *

def preprocess_data(session_ID, adata,
                    min_cells=2, min_genes=200, max_genes=10000,
                    target_sum=1e6, flavor="cell_ranger", 
                    n_top_genes=2000):

    adata = generate_adata_from_10X(session_ID)

    # do preprocessing
    print("[STATUS] performing QC and normalizing data")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    #adata.raw = adata


    # filter down to top 2000 highly-variable genes, using cell_ranger method
    # TODO: parametrize this in the webapp
    print("[STATUS] selecting highly variable genes")
    sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top_genes)
    #adata = adata[:, adata.var['highly_variable']]

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

def do_pseudotime(session_ID, adata, starter_cell_ID=None):
    a = adata.to_df()

    print("[STATUS] computing pseudotime ...")
    pca_projections, var_r = palantir.utils.run_pca(a)
    #adata.uns['palantir_pca_results'] = {}
    #adata.uns['palantir_pca_results']['pca_projections'] = pca_projections
    #adata.uns['palantir_pca_results']['variance_ratio']  = var_r

    dm_res = palantir.utils.run_diffusion_maps(pca_projections, knn=20)
    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    #adata.uns['palantir_diff_maps'] = dm_res
    #adata.uns['palantir_ms_data'] = ms_data

    #d.tsne = d.palantir.utils.run_tsne(d.ms_data, perplexity=150)
    # We're just going to use the UMAP projection instead
    tsne = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs.index)
    tsne.columns = ["x", "y"]
    #adata.uns['palantir_tsne'] = tsne

    imp_df = palantir.utils.run_magic_imputation(a, dm_res, n_steps=1)
    #adata.uns['palantir_imp_df'] = imp_df

    
    #start_cell = str(adata.obs.index[0])
    if (starter_cell_ID is None):
        start_cell = "TTGTTTGCAATTTCCT"
        print("[ERROR] no starter cell provided; using " + start_cell +
              " as a default (assuming original type-II data)")
    else:
        start_cell = starter_cell_ID
    #start_cell = "TTGTTTGCAATTTCCT-1"
    #start_cell = "ATGGAGGCAGCTAACT-1" #fails
    print("[DEBUG] start cell is: " + str(start_cell))
    pr_res = palantir.core.run_palantir(ms_data, start_cell, 
                                        terminal_states=None, knn=20, 
                                        num_waypoints=500, 
                                        use_early_cell_as_start=True, 
                                        scale_components=False)
    adata.obs["pseudotime"] = pr_res.pseudotime[tsne.index]
    adata.obs["differentiation_potential"] = pr_res.entropy[tsne.index]
    #adata.uns["pr_res"] = pr_res

    genes = adata.var.index.tolist()
    #genes = ["CycE", "dap", "Hey", "nSyb"]
    #print("[DEBUG] genes: " + str(genes))
    #print("[DEBUG] d.imp_df.loc[:, genes]: " + str(d.imp_df.loc[:, genes]))
    print("[STATUS] computing all gene trends (this will take a while)")
    gene_trends = palantir.presults.compute_gene_trends(pr_res, 
                                                        imp_df.loc[:, genes])

    #d.adata.uns["gene_trends"] = gene_trends
    cache_adata(session_ID, adata)
    cache_gene_trends(session_ID, gene_trends)
    return adata

def do_pseudotime_gene_trends(session_ID, d, genes):
    d = cache_pseudotime(session_ID)
    gene_trends = d.palantir.presults.compute_gene_trends(d.adata.uns["pr_res"], 
                                                          d.imp_df.loc[:, genes])
    cache_gene_trends(session_ID, gene_trends)
    return gene_trends