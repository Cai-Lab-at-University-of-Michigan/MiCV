import scanpy as sc

def identify_marker_genes(adata, obs_column, groups_to_rank, method):
    if ("all" in groups_to_rank):
        sc.tl.rank_genes_groups(adata, groupby=obs_column, method=method)
    else:
        sc.tl.rank_genes_groups(adata, groupby=obs_column, 
                                groups=groups_to_rank, method=method)
    return adata