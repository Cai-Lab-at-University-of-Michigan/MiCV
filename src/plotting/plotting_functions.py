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
from plotting.multi_color_scale import MultiColorScale
from plotting.discrete_color_scales import *
from plotting.plotting_parameters import scale

min_opacity = 0.2
max_opacity = 0.9

def plot_UMAP(adata, clustering_plot_type, selected_cell_intersection=[], n_dim=2):
    print("[DEBUG] generating new UMAP plot")

    # validate that there is a 3D projection available if that was requested
    if (n_dim == 3):
        if not ("X_umap_3D" in adata.obsm):
            n_dim = 2
    
    traces = []
    for i,val in enumerate(sorted(adata.obs[clustering_plot_type].unique())):
        a = adata[adata.obs[clustering_plot_type] == val]
        s = []
        if (selected_cell_intersection in [None, []]):
        	s = list(range(0, len((a.obs).index)))
        else:
	        for c in selected_cell_intersection:
	            if (c in a.obs["cell_numeric_index"]):
	                s.append((a.obs).index.get_loc(c))
        if (n_dim == 2):
            traces.append(
                go.Scattergl(
                    x=a.obsm["X_umap"][:,0],
                    y=a.obsm["X_umap"][:,1],
                    text="Cell ID: " + a.obs["cell_ID"],
                    mode='markers',
                    selectedpoints=s,
                    marker={
                        'size': 10,
                        'line': {'width': 1, 'color': 'grey'},
                        "color": discrete_colors_3[i%len(discrete_colors_3)]
                    },
                    unselected={
                        "marker": {"opacity": min_opacity,
                        }
                    },
                    selected={
                        "marker": {"opacity": max_opacity,
                        }
                    },
                    name=("Cluster " + str(val))
                )
            )
        elif (n_dim == 3):
            traces.append(
                go.Scatter3d(
                    x=a.obsm["X_umap_3D"][:,0],
                    y=a.obsm["X_umap_3D"][:,1],
                    z=a.obsm["X_umap_3D"][:,2],
                    text="Cell ID: " + a.obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': 2.5,
                        'line': {'width': 1, 'color': 'grey'},
                        "color": discrete_colors_3[i%len(discrete_colors_3)]
                    },
                    name=("Cluster " + str(val))
                )
            )
    if (n_dim == 2):
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
    elif (n_dim == 3):
        return {
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                zaxis={"title": "UMAP 3"},
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
    elif ("pseudotime_branch_" in pt_plot_type):
        colorbar_label = "branch " + str(pt_plot_type[-1]) + " prob."
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

def plot_expression_UMAP(adata, selected_genes, multi="standard", n_dim=2):
    
    # validate that there is a 3D projection available if that was requested
    if (n_dim == 3):
        if not ("X_umap_3D" in adata.obsm):
            n_dim = 2

    if (multi == "standard"):
        colorscale = "viridis"
        selected_gene = selected_genes
        traces = []
        if (n_dim == 2):
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
                        "colorscale": colorscale,
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
        elif (n_dim == 3):
            traces.append(
                go.Scatter3d(
                    x=adata.obsm["X_umap_3D"][:,0],
                    y=adata.obsm["X_umap_3D"][:,1],
                    z=adata.obsm["X_umap_3D"][:,2],
                    text="Cell ID: " + adata.obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': 2.5,
                        'line': {'width': 1, 'color': 'grey'},
                        "color": adata.obs_vector(selected_gene),
                        "colorscale": colorscale,
                        "cmin": 0,
                        "cmax": np.max(adata.obs_vector(selected_gene)),
                        "colorbar": dict(
                            title=str(selected_gene)
                        ),
                    },
                )
            )
    
    else:
        if (len(selected_genes) > 3):
            selected_genes = selected_genes[0:3]
        color_values = get_mixed_expression_value(*[adata.obs_vector(selected_genes[g]) for g in range(0, len(selected_genes))])
        traces = []
        if (n_dim == 2):
            traces.append(
                go.Scattergl(
                    x=adata.obsm["X_umap"][:,0],
                    y=adata.obsm["X_umap"][:,1],
                    text="Cell ID: " + adata.obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': 10,
                        'line': {'width': 1, 'color': 'grey'},
                        "color": color_values,
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
        elif (n_dim == 3):
            traces.append(
                go.Scatter3d(
                    x=adata.obsm["X_umap_3D"][:,0],
                    y=adata.obsm["X_umap_3D"][:,1],
                    z=adata.obsm["X_umap_3D"][:,2],
                    text="Cell ID: " + adata.obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': 2.5,
                        'line': {'width': 1, 'color': 'grey'},
                        "color": color_values,
                    },
                )
            )    

    if (n_dim == 2):
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
    elif (n_dim == 3):
        return {
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                zaxis={"title": "UMAP 3"},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                width=4 * scale,
                height=3 * scale
            )
        }

def plot_expression_trend(gene_trends, selected_genes, selected_branch, 
                          relative="absolute"):

    traces = []
    trends = gene_trends[selected_branch]["trends"]
    stds = gene_trends[selected_branch]["std"] * 25

    colors = pd.Series(sns.color_palette('Set2', len(selected_genes)).as_hex(), 
                       index=selected_genes)
    for i in selected_genes:
        if not (i in trends.index):
            print("[DEBUG] gene " + str(i)  + " not in gene trends; skipping")
            continue
        if (relative == "relative"):
            trend = trends.loc[i,:] / np.max(trends.loc[i,:])
            std   = stds.loc[i,:] / np.max(stds.loc[i,:])
        else:
            trend = trends.loc[i,:]
            std   = stds.loc[i,:]
        traces.append(
            go.Scattergl(
                x=trends.columns,
                y=trend,
                text=str(i),
                mode="lines+markers",
                opacity=0.7,
                marker={
                    'size': std,
                    'line': {'width': 2, 'color': colors[i]},
                    "color": colors[i],
                    "opacity": 0.25
                },
                name=(str(i))
            )
        )
    
    if (traces in [[], None]):
        print("[DEBUG] no traces added to expression trends plot")
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
            width=4 * scale,
            height=3 * scale
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
            width=4 * scale,
            height=3 * scale
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

def get_mixed_expression_value(e0, e1=None, e2=None):
    m = 1
    n_levels = 16

    def rescale_expression(a, m):
        a = np.asarray(a)
        a = m * (a - np.min(a))/np.ptp(a)
        return a

    e0 = rescale_expression(e0, m)

    if (e1 is None):
        e1 = np.zeros_like(e0)
    else:
        e1 = rescale_expression(e1, m)
    if (e2 is None):
        e2 = np.zeros_like(e0)
    else:
        e2 = rescale_expression(e2, m)

    colors = np.column_stack((e0,e1,e2))
    print("[DEBUG] colors: " + str(colors[0:5]))

        
    bins = np.linspace(0, m, n_levels)
    scale = cache_multicolor_scale()

    inds = np.digitize(colors, bins, right=True)
    for i in range(0, len(colors)):
        for j in range(0, len(colors[i])):
            colors[i,j] = bins[min(inds[i,j], n_levels - 1)]

    ret = [scale.calculate_hex_color(c) for c in colors]
    return ret