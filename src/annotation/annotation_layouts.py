import dash_html_components as html
import dash_bootstrap_components as dbc

from . import annotation_components as cc #custom components

def annotation_layout():
	m = dbc.Tab(label='Annotation', children=[
			    html.Div(children=[
				    dbc.Button("Load processed data", 
				    			id="load_analysis_button"),
		    		html.P("Press this button to get started by loading the processed data")
				]),  

			    html.Div(children=[
				    dbc.Button("Define new cluster", 
				    			id="define_cluster_button"),
				    html.P("Add currently selected cells from the clustering plot to a new cluster in the user-cluster plot group you have currently selected", id="define_cluster_text")
			    ],
			    style={'marginBottom': 20, 'marginTop': 20}),	
			    
			    # Plots
				dbc.Row(children=[
					dbc.Col(children=[
					   	# clustering plot
					    html.H3(children="Clustering plot"),
					    html.P("Show either leiden (automatic) cluster assignments or user-defined cluster assignments"),
					    cc.clustering_dropdown(),
					    cc.plot_clustering_UMAP(),
					    html.H5(children=[cc.clustering_UMAP_count()]),
					    ], width=6,
					),
				    dbc.Col(children=[
					   	# pseudotime plot
					    html.H3(children="Pseudotime plot"),
					    html.P("Show either the pseudotime assignment or calculated differentiation potential of cells, based on the palantir markov-chain-based algorithm"),
					    cc.pseudotime_dropdown(),
					    cc.plot_pseudotime_UMAP(),
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
					    cc.plot_expression_UMAP(),
					    html.H5(children=[cc.gene_UMAP_count()]),
					    html.H3(children="Gene information"),
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
					    cc.plot_gene_pseudotime(),
					    cc.plot_gene_violin(),
					    html.H5(children=[cc.gene_violin_count()]),
					    ], width=6,
					)
			    ], 
			    )
			]) ## end annotation tab
	return m