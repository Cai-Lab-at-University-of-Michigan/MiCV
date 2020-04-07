import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

from . import annotation_components as cc #custom components

def annotation_layout():
	m = dbc.Tab(label='Annotation', children=[ 
		    	dbc.Col(children=[
		    		html.H3(children="Data exploration and expert annotation"),
					html.P(
						'''
						This is the main interface for the MiCV tool.
						Use the plots below to explore automated clustering,
						gene expression, gene literature references, and
						pseudotime trajectories/expression dynamics (if
						calculated previously). Then, use the selection tools 
						on any and all of the plots to select cells of interest
						based on your expert analysis, and define a new manually-
						annotated cluster by selecting one of the user_n cluster
						options in the top-left clustering plot dropdown menu and
						pressing the define new cluster button below. 
						''', id="define_cluster_text", style={'marginBottom': 10, 
															  'marginTop': 10}),
			    	dbc.Button("Define new cluster", 
			    			id="define_cluster_button"),
		    	], width=6),
			    
			    # Plots
				dbc.Row(children=[
					dbc.Col(children=[
					   	# clustering plot
					    html.H3(children="Clustering plot"),
					    html.P("Show either leiden (automatic) cluster assignments or user-defined cluster assignments"),
					    cc.clustering_dropdown(),
					    dcc.Loading(children=[cc.plot_clustering_UMAP()]),
					    html.H5(children=[cc.clustering_UMAP_count()]),
					    ], width=6,
					),
				    dbc.Col(children=[
					   	# pseudotime plot
					    html.H3(children="Pseudotime plot"),
					    html.P("Show either the pseudotime assignment or calculated differentiation potential of cells, based on the palantir markov-chain-based algorithm"),
					    cc.pseudotime_dropdown(),
					    dcc.Loading(children=[cc.plot_pseudotime_UMAP()]),
					    html.H5(children=[cc.pseudotime_UMAP_count()]),
					    cc.total_cell_count()
					    ], width=6,
					)
			    ], 
			    ),

			    dbc.Row(children=[
					dbc.Col(children=[
					   	# expression plot
					    html.H3(children="Gene expression (projection)"),
					    html.P("Visualize expression of a single highly-variable gene across single cells, and pull up data from Flybase (Oct, 2019) on the selected gene"),
					    dbc.Row(children=[
					    	dbc.Col(children=[
					    		cc.single_gene_dropdown(),
					    		cc.mixed_gene_dropdown()
					    	], width=9),
					    	dbc.Col(children=[
					    		cc.single_gene_expression_radio(),
					    		cc.n_dims_proj_expression_radio()
					    	], width=3),
						]),
					    dbc.Row(children=[
					    	dbc.Col(children=[
					    		dcc.Loading(children=[cc.plot_expression_UMAP()]),
					    	], width=9),
					    	dbc.Col(children=[
					    		html.Img(src="assets/color_triangle_legend.png", 
					    				 hidden=True,
					    				 id="mixed_expression_legend_image",
					    				 style={"display": "flex",
												"justify-content": "center",
												"align-items": "center",
												"width": "90%"})
					    	], width=3)
					    ]),
					    html.H5(children=[cc.gene_UMAP_count()]),
					    html.H3(children="Gene information (fly only)"),
					    cc.gene_data_table()
					    ], width=6,
					),
				    dbc.Col(children=[
					   	# pseudotime gene expression plot
					    html.H3(children="Gene expression (pseudotime/bulk)"),
					    html.P("Visualize the expression of multiple highly-variable genes against pseudotime or across all cells from the entire dataset"),
					    cc.multi_gene_dropdown(),
					    cc.pseudotime_gene_relative_radio(),
					    cc.pseudotime_gene_branch_dropdown(),
					    dcc.Loading(children=[cc.plot_gene_pseudotime()]),
					    dcc.Loading(children=[cc.plot_gene_violin()]),
					    html.H5(children=[cc.gene_violin_count()]),
					    ], width=6,
					)
			    ], 
			    )
			]) ## end annotation tab
	return m