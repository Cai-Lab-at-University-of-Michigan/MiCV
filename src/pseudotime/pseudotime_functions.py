import pandas as pd
import numpy as np
from pygam import LinearGAM, ExpectileGAM, s, f
import palantir

from helper_functions import *
from status.status_functions import *

def do_pseudotime(session_ID, adata, starter_cell_ID=None):
    n_steps = 6
    a = adata.to_df()
    cache_progress(session_ID, progress=int(1/n_steps * 100))
    print("[STATUS] computing pseudotime ...")

    pca_projections, var_r = palantir.utils.run_pca(a)
    cache_progress(session_ID, progress=int(2/n_steps * 100))
    cache_history(session_ID, history=("Performed PCA for pseudotime"))

    dm_res = palantir.utils.run_diffusion_maps(pca_projections, knn=20)
    cache_progress(session_ID, progress=int(3/n_steps * 100))

    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    cache_progress(session_ID, progress=int(4/n_steps * 100))
    cache_history(session_ID, history=("Generated diffusion maps and"
                                     + " multiscale space for pseudotime"))

    imp_df = palantir.utils.run_magic_imputation(a, dm_res, n_steps=1)
    cache_progress(session_ID, progress=int(5/n_steps * 100))
    cache_history(session_ID, history=("Imputed data for pseudotime"))
    
    start_cell = starter_cell_ID
    # TODO: need to catch ValueError here in try-catch block
    pr_res = palantir.core.run_palantir(ms_data, start_cell, 
                                        terminal_states=None,
                                        num_waypoints=500, 
                                        scale_components=True, n_jobs=1)
    cache_progress(session_ID, progress=int(5/n_steps * 100))
    cache_history(session_ID, history=("Calculated pseudotime trajectories with "
                                     + "starter cell ID: " + str(starter_cell_ID)))

    adata.obs["pseudotime"] = pr_res.pseudotime[adata.obs.index]
    adata.obs["differentiation_potential"] = pr_res.entropy[adata.obs.index]
    
    # drop old branch probabilities if they exist, then add new ones
    adata.obs.drop(list(adata.obs.filter(regex='pseudotime_branch_')), 
                   axis=1, inplace=True)
    for i, branch in enumerate(pr_res.branch_probs.columns):
        adata.obs["pseudotime_branch_" + str(i)] = pr_res.branch_probs.loc[adata.obs.index, branch]
        #adata.obs["in_pseudotime_branch_" + str(i)] = adata.obs["pseudotime_branch_" + str(i)] >= 0.7

    cache_adata(session_ID, adata)
    cache_pseudotime_results(session_ID, pr_res)
    cache_imputed_df(session_ID, imp_df)
    cache_progress(session_ID, progress=int(6/n_steps * 100))
    return adata

def calculate_gene_trends(session_ID, list_of_genes, branch_ID):
    n_steps = 2 + len(list_of_genes)
    
    pr_res = cache_pseudotime_results(session_ID)
    cache_progress(session_ID, progress=int(1/n_steps * 100))

    imp_df = cache_imputed_df(session_ID)
    obs = cache_adata(session_ID, group="obs")
    cache_progress(session_ID, progress=int(2/n_steps * 100))

    if (branch_ID == -1):
        branch = "all branches"
        cells_in_branch = obs.index
    else:
        branch = list(pr_res.branch_probs.columns)[branch_ID]
        cells_in_branch = obs[obs["pseudotime_branch_" + str(branch_ID)] > 0.7].index
    print("[DEBUG] branch: " + str(branch))

    '''
    gene_trends = palantir.presults.compute_gene_trends(pr_res, 
                                                        imp_df.loc[:, genes],
                                                        lineages = [branch],
                                                        n_jobs=1)
    '''
    X_train = pr_res.pseudotime[cells_in_branch].to_numpy()
    X_train = np.reshape(X_train, (len(X_train), 1))

    X_plot = np.linspace(np.min(X_train), np.max(X_train), 50)
    
    gene_trends = pd.DataFrame()
    gene_trends["pseudotime"] = X_plot
    
    step_number = 3
    for gene in list_of_genes:
        Y_train = imp_df.loc[cells_in_branch, gene].to_numpy()
        
        #gam = LinearGAM(n_splines=25)
        #gam.gridsearch(X_train, Y_train, progress=False)
        gam = ExpectileGAM(expectile=0.5).gridsearch(X_train, Y_train)
        lam = gam.lam
        gam_upper = ExpectileGAM(expectile=0.95, lam=lam).fit(X_train, Y_train)
        gam_lower = ExpectileGAM(expectile=0.05, lam=lam).fit(X_train, Y_train)
        gene_trends[gene] = gam.predict(X_plot)
        
        #ci = gam.prediction_intervals(X_plot, width=.95)
        gene_trends[gene + "_ci_upper"] = gam_upper.predict(X_plot)
        gene_trends[gene + "_ci_lower"] = gam_lower.predict(X_plot)
        cache_progress(session_ID, progress=int(step_number/n_steps * 100))
        step_number += 1
    return gene_trends