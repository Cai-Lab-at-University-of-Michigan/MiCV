import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import anndata as ad
import seaborn as sns

from scipy.stats import gaussian_kde

import base64

from helper_functions import *
from plotting.multi_color_scale import MultiColorScale
from plotting.discrete_color_scales import *
from plotting.plotting_parameters import *


def plot_UMAP(session_ID, clustering_plot_type, selected_cell_intersection=[], n_dim=2):
    print("[DEBUG] generating new UMAP plot")

    if (adata_cache_exists(session_ID) is False):
        print("[ERROR] cache for " + str(session_ID) + "does not exist")
        return dash.no_update
    if ((adata_cache_group_exists(session_ID, "obs") is False)
    or  (adata_cache_group_exists(session_ID, "obsm") is False)):
        print("[ERROR] obs/obsm for " + str(session_ID) + "does not exist")
        return dash.no_update

    obs  = cache_adata(session_ID, group="obs")
    obsm = cache_adata(session_ID, group="obsm")

    # validate that there is a 3D projection available if that was requested
    if (("X_umap_3D" in obsm.keys()) and (n_dim == 3)):
        coords = pd.DataFrame(obsm["X_umap_3D"], index=obs.index)
    else:
        n_dim = 2
        coords = pd.DataFrame(obsm["X_umap"], index=obs.index)

    
    traces = []
    for i,val in enumerate(sorted(obs[clustering_plot_type].unique())):
        a = obs[obs[clustering_plot_type] == val]
        b = coords[obs[clustering_plot_type] == val]
        s = []
        if (selected_cell_intersection in [None, []]):
        	s = list(range(0, len(a.index)))
        else:
	        for c in selected_cell_intersection:
	            if (c in a["cell_numeric_index"]):
	                s.append(a.index.get_loc(c))
        if (n_dim == 2):
            traces.append(
                go.Scattergl(
                    x=b[0],
                    y=b[1],
                    text="Cell ID: " + a["cell_ID"],
                    mode='markers',
                    selectedpoints=s,
                    marker={
                        'size': point_size_2d,
                        'line': {'width': point_line_width_2d, 'color': 'grey'},
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
                    x=b[0],
                    y=b[1],
                    z=b[2],
                    text="Cell ID: " + a["cell_ID"],
                    mode='markers',
                    marker={
                        'size': point_size_3d,
                        'line': {'width': point_line_width_3d, 'color': 'grey'},
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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        }
    elif (n_dim == 3):
        return {
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                zaxis={"title": "UMAP 3"},
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        }

def plot_pseudotime_UMAP(session_ID, pt_plot_type, n_dim=2):
    if (adata_cache_exists(session_ID) is False):
        return dash.no_update
    if ((adata_cache_group_exists(session_ID, "obs") is False)
    or  (adata_cache_group_exists(session_ID, "obsm") is False)):
        print("[ERROR] obs/obsm for " + str(session_ID) + "does not exist")
        return dash.no_update

    obs  = cache_adata(session_ID, group="obs")
    obsm = cache_adata(session_ID, group="obsm")

    if not (pt_plot_type in obs):
        return dash.no_update
    
    # validate that there is a 3D projection available if that was requested
    if (("X_umap_3D" in obsm.keys()) and (n_dim == 3)):
        coords = pd.DataFrame(obsm["X_umap_3D"], index=obs.index)
    else:
        n_dim = 2
        coords = pd.DataFrame(obsm["X_umap"], index=obs.index)

    if (pt_plot_type == "pseudotime"):
    	colorbar_label = "pseudotime"
    elif (pt_plot_type == "differentiation_potential"):
    	colorbar_label = "diff. pot."
    elif ("pseudotime_branch_" in pt_plot_type):
        colorbar_label = "branch " + str(pt_plot_type[-1]) + " prob."

    traces = []
    traces.append(
        go.Scattergl(
            x=coords[0],
            y=coords[1],
            text="Cell ID: " + obs["cell_ID"],
            mode='markers',
            marker={
                'size': point_size_2d,
                'line': {'width': point_line_width_2d, 'color': 'grey'},
                "color": obs[str(pt_plot_type)],
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
            margin=margin,
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 250},
            autosize=True
            #width=4 * scale,
            #height=3 * scale
        )
    }

def plot_expression_UMAP(session_ID, selected_genes, multi="standard", n_dim=2):
    
    if (adata_cache_exists(session_ID) is False):
        return dash.no_update

    adata  = cache_adata(session_ID)
    obsm   = adata.obsm
    obs    = adata.obs

    # validate that there is a 3D projection available if that was requested
    if (("X_umap_3D" in obsm.keys()) and (n_dim == 3)):
        coords = pd.DataFrame(obsm["X_umap_3D"], index=obs.index)
    else:
        n_dim = 2
        coords = pd.DataFrame(obsm["X_umap"], index=obs.index)

    
    if (multi == "standard"):
        colorscale = "viridis"
        selected_gene = selected_genes
        traces = []
        if (n_dim == 2):
            traces.append(
                go.Scattergl(
                    x=coords[0],
                    y=coords[1],
                    text="Cell ID: " + obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': point_size_2d,
                        'line': {'width': point_line_width_2d, 'color': 'grey'},
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
                    x=coords[0],
                    y=coords[1],
                    z=coords[2],
                    text="Cell ID: " + obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': point_size_3d,
                        'line': {'width': point_line_width_3d, 'color': 'grey'},
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
                    x=coords[0],
                    y=coords[1],
                    text="Cell ID: " + obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': point_size_2d,
                        'line': {'width': point_line_width_2d, 'color': 'grey'},
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
                    x=coords[0],
                    y=coords[1],
                    z=coords[2],
                    text="Cell ID: " + obs["cell_ID"],
                    mode='markers',
                    marker={
                        'size': point_size_3d,
                        'line': {'width': point_line_width_3d, 'color': 'grey'},
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
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        }
    elif (n_dim == 3):
        return {
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                zaxis={"title": "UMAP 3"},
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        }

def plot_expression_trend(gene_trends, selected_genes, selected_branch, 
                          relative="absolute"):

    traces = []
    trends = gene_trends

    colors = pd.Series(sns.color_palette('Set2', len(selected_genes)), 
                       index=selected_genes)
    fill_opacity = 0.1
    '''
    fill_colors = list(sns.color_palette('Set2', len(selected_genes)))
    print(fill_colors)
    for i, c in enumerate(fill_colors):
        fill_colors[i] = tuple(list(fill_colors[i]).append(fill_opacity))
    fill_colors = pd.Series(fill_colors, index=selected_genes)
    '''
    for i in selected_genes:
        if not (i in trends.columns):
            print("[DEBUG] gene " + str(i)  + " not in gene trends; skipping")
            continue
        if (relative == "relative"):
            trend = trends[i] / np.max(trends[i])
            ci_upper = trends[i+"_ci_upper"] / np.max(trends[i])
            ci_lower = trends[i+"_ci_lower"] / np.max(trends[i])
        else:
            trend = trends[i]
            ci_upper = trends[i+"_ci_upper"]
            ci_lower = trends[i+"_ci_lower"]
        X = trends["pseudotime"]
        traces.append(
            go.Scatter(
                x=X,
                y=ci_upper,
                showlegend=False,
                mode="lines",
                line_color=to_rgba_string(colors[i], fill_opacity),
                fill=None,
                name=str(i)
            )
        )
        traces.append(
            go.Scatter(
                x=X,
                y=ci_lower,
                showlegend=False,
                fill='tonexty',
                mode="lines",
                line_color=to_rgba_string(colors[i], fill_opacity),
                fillcolor=to_rgba_string(colors[i], fill_opacity),
                name=str(i)
            )
        )
        traces.append(
            go.Scatter(
                x=X,
                y=trend,
                text=str(i),
                mode="markers+lines",
                opacity=1,
                name=(str(i)),
                marker={
                    'size': point_size_pt_trend,
                },
                line_color=to_rgba_string(colors[i])
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
            margin=margin,
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
            autosize=True
        )
    }


def plot_expression_violin(session_ID, selected_genes, show_points = "all"):
    if (adata_cache_exists(session_ID) is False):
        return dash.no_update
    
    adata  = cache_adata(session_ID)
    var    = adata.var
    obs    = adata.obs

    traces = []
    
    x_pos  = 1
    n_traces = len(selected_genes)
    
    for i in selected_genes:
        if not ((i in var.index) or (i in obs)):
            print("[DEBUG] gene " + str(i) + " not in var index; skipping")
            continue
        if (show_points == False):
            traces.append(
                go.Violin(
                    y=adata.obs_vector(i),
                    text="Cell ID: " + obs["cell_ID"],
                    opacity=0.7,
                    box_visible=True,
                    meanline_visible=True,
                    points=False,
                    name=(str(i))
                )
            )
            x_pos += 1
        
        elif (show_points == "all"):
            #kernel = gaussian_kde(adata.obs_vector(i))
            jittered_x = x_pos + 0.1 * np.random.standard_normal(len(adata.obs_vector(i)))

            traces.append(
                go.Scattergl(
                    x=jittered_x,
                    y=adata.obs_vector(i),
                    text="Cell ID: " + obs["cell_ID"],
                    mode="markers",
                    opacity=0.7,
                    marker={
                        'size': point_size_2d,
                    },
                    name=(str(i)),
                )
            )
            x_pos += 1

    if (traces in [[], None]):
        print("[DEBUG] no traces added to violin plot")
        return dash.no_update

    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Gene/factor"},
            yaxis={"title": "Expression"},
            margin=margin,
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
            autosize=True
            #width=4 * scale,
            #height=3 * scale
        )
    }

def plot_marker_genes(session_ID, adata, obs_column, groups_to_rank):
    sc.settings.figdir = save_analysis_path + str(session_ID) + "/"
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
    encoded_image = base64.b64encode(open(save_analysis_path + str(session_ID) + "/" + image_filename, 'rb').read())
    return html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()), style={"width": "95%"})

def get_mixed_expression_value(e0, e1=None, e2=None, session_ID=None):
    m = 1
    
    def rescale_expression(a, m):
        a = np.asarray(a)
        a = m * (a - np.min(a))/np.ptp(a)
        return a

    def expression_to_rgb(e):
        ret = np.column_stack((e[:,1]+e[:,2],
                               e[:,0]+e[:,2],
                               e[:,0]+e[:,1]))
        return ret

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
    ret = expression_to_rgb(colors)

    return ret
