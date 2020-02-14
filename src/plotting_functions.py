import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import anndata as ad
import seaborn as sns

import base64

from helper_functions import *
from markergenes import *

from plotting.plotting_parameters import scale, min_opacity, max_opacity

def plot_UMAP(adata, clustering_plot_type, selected_cell_intersection=[]):
    print("[DEBUG] generating new UMAP plot")

    traces = []
    for i in sorted(adata.obs[clustering_plot_type].unique()):
        a = adata[adata.obs[clustering_plot_type] == i]
        s = []
        if (selected_cell_intersection in [None, []]):
        	s = list(range(0, len((a.obs).index)))
        else:
	        for c in selected_cell_intersection:
	            if (c in a.obs["cell_numeric_index"]):
	                s.append((a.obs).index.get_loc(c))
        traces.append(
            go.Scattergl(
                x=a.obsm["X_umap"][:,0],
                y=a.obsm["X_umap"][:,1],
                text="Cell ID: " + a.obs["cell_ID"],
                mode='markers',
                selectedpoints=s,
                marker={
                    'size': 10,
                    'line': {'width': 1, 'color': 'grey'}
                },
                unselected={
                    "marker": {"opacity": min_opacity,
                    }
                },
                selected={
                    "marker": {"opacity": max_opacity,
                    }
                },
                name=("Cluster " + str(i))
            )
        )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "UMAP 1"},
            yaxis={"title": "UMAP 2"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
            width=4 * scale,
            height=3 * scale
        )
    }

def plot_pseudotime_UMAP(adata, pt_plot_type):
    if (pt_plot_type == "pseudotime"):
    	colorbar_label = "pseudotime"
    elif (pt_plot_type == "differentiation_potential"):
    	colorbar_label = "diff. pot."
    traces = []
    traces.append(
        go.Scattergl(
            x=adata.obsm["X_umap"][:,0],
            y=adata.obsm["X_umap"][:,1],
            text="Cell ID: " + adata.obs["cell_ID"],
            mode='markers',
            marker={
                'size': 10,
                'line': {'width': 1, 'color': 'grey'},
                "color": adata.obs[str(pt_plot_type)],
                "colorscale": "plasma",
                "cmin": 0,
                "cmax": 1,
                "colorbar": dict(
                    title=colorbar_label
                ),
            },
            unselected={
                "marker": {"opacity": min_opacity,
                }
            },
            selected={
                "marker": {"opacity": max_opacity,
                }
            },
        )
    )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "UMAP 1"},
            yaxis={"title": "UMAP 2"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
            width=4 * scale,
            height=3 * scale
        )
    }

def plot_expression_UMAP(adata, selected_gene):
    traces = []
    traces.append(
        go.Scattergl(
            x=adata.obsm["X_umap"][:,0],
            y=adata.obsm["X_umap"][:,1],
            text="Cell ID: " + adata.obs["cell_ID"],
            mode='markers',
            marker={
                'size': 10,
                'line': {'width': 1, 'color': 'grey'},
                "color": adata.obs_vector(selected_gene),
                "colorscale": "viridis",
                "cmin": 0,
                "cmax": np.max(adata.obs_vector(selected_gene)),
                "colorbar": dict(
                    title=str(selected_gene)
                ),
            },
            unselected={
                "marker": {"opacity": min_opacity,
                }
            },
            selected={
                "marker": {"opacity": max_opacity,
                }
            },
        )
    )
    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "UMAP 1"},
            yaxis={"title": "UMAP 2"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
            width=4 * scale,
            height=3 * scale
        )
    }

def plot_expression_trend(gene_trends, selected_branch, selected_genes):
    traces = []
    trends = gene_trends[selected_branch]["trends"]
    stds = gene_trends[selected_branch]["std"] * 25
    colors = pd.Series(sns.color_palette('Set2', len(selected_genes)).as_hex(), 
                       index=selected_genes)
    for i in selected_genes:
        if not (i in trends.index):
            print("[DEBUG] gene " + str(i)  + " not in gene trends; skipping")
            continue
        traces.append(
            go.Scattergl(
                x=trends.columns,
                y=trends.loc[i, :],
                text=str(i),
                mode="lines+markers",
                opacity=0.7,
                marker={
                    'size': stds.loc[i, :],
                    'line': {'width': 2, 'color': colors[i]},
                    "color": colors[i],
                    "opacity": 0.25
                },
                name=(str(i))
            )
        )
    
    if (traces in [[], None]):
        print("[DEBUG] no traces added to violin plot")
        return dash.no_update

    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Pseudotime"},
            yaxis={"title": "Expression"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
        )
    }

def plot_expression_violin(adata, selected_genes):
    traces = []
    for i in selected_genes:
        if not ((i in adata.var.index) or (i in adata.obs)):
            print("[DEBUG] gene " + str(i) + " not in var index; skipping")
            continue
        d = adata.obs_vector(i)
        traces.append(
            go.Violin(
                y=d,
                text="Cell ID: " + adata.obs["cell_ID"],
                opacity=0.7,
                box_visible=True,
                meanline_visible=True,
                points="all",
                name=(str(i))
            )
        )

    if (traces in [[], None]):
        print("[DEBUG] no traces added to violin plot")
        return dash.no_update

    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Gene/factor"},
            yaxis={"title": "Expression"},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
        )
    }

def plot_marker_genes(adata, obs_column, groups_to_rank):
    sc.settings.figdir = save_analysis_path
    print("[STATUS] identifying marker genes")
    image_filename = "dotplot.png"
    
    print("[STATUS] plotting marker genes")
    if ("all" in groups_to_rank):
        ax = sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups", 
                                         dendrogram=False, groupby=obs_column, show=False,
                                         save=".png")
    else:
        ax = sc.pl.rank_genes_groups_dotplot(adata[adata.obs[obs_column].isin(groups_to_rank),:], 
                                             groups=groups_to_rank, key="rank_genes_groups", 
                                             dendrogram=False, groupby=obs_column, show=False,
                                             save=".png")
    encoded_image = base64.b64encode(open(save_analysis_path + image_filename, 'rb').read())
    return html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()))