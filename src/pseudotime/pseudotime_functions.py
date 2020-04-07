import pandas as pd
import palantir

from helper_functions import *

def do_pseudotime(session_ID, adata, starter_cell_ID=None):
    a = adata.to_df()

    print("[STATUS] computing pseudotime ...")
    pca_projections, var_r = palantir.utils.run_pca(a)

    dm_res = palantir.utils.run_diffusion_maps(pca_projections, knn=20)
    ms_data = palantir.utils.determine_multiscale_space(dm_res)

    imp_df = palantir.utils.run_magic_imputation(a, dm_res, n_steps=1)
    
    start_cell = starter_cell_ID
    
    pr_res = palantir.core.run_palantir(ms_data, start_cell, 
                                        terminal_states=None,
                                        num_waypoints=500, 
                                        scale_components=True, n_jobs=1)
    adata.obs["pseudotime"] = pr_res.pseudotime[adata.obs.index]
    adata.obs["differentiation_potential"] = pr_res.entropy[adata.obs.index]
    
    # drop old branch probabilities if they exist, then add new ones
    adata.obs.drop(list(adata.obs.filter(regex='pseudotime_branch_')), 
                   axis=1, inplace=True)
    for i, branch in enumerate(pr_res.branch_probs.columns):
        adata.obs["pseudotime_branch_" + str(i)] = pr_res.branch_probs.loc[adata.obs.index, branch]

    cache_adata(session_ID, adata)
    cache_pseudotime_results(session_ID, pr_res)
    cache_imputed_df(session_ID, imp_df)
    return adata

def calculate_gene_trends(session_ID, list_of_genes, branch_ID):
    pr_res = cache_pseudotime_results(session_ID)
    imp_df = cache_imputed_df(session_ID)

    branch = list(pr_res.branch_probs.columns)[branch_ID]
    print("[DEBUG] branch: " + str(branch))

    genes = list_of_genes
    gene_trends = palantir.presults.compute_gene_trends(pr_res, 
                                                        imp_df.loc[:, genes],
                                                        lineages = [branch],
                                                        n_jobs=1)
    return gene_trends