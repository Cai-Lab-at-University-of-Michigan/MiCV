import dash_core_components as dcc
import dash_html_components as html
import custom_components as cc
import uuid

def main_layout():
	session_id = str(uuid.uuid4())
	#session_id = "1234" # don't use in production
	return html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    
	    # title matter
	    html.H1(children="MiCV"),
	    html.H3(children="A Multi-Informatic Cellular Visualization tool"),

	    # clustering parameters
	    html.Div(id="analysis_page", children=[
	    	html.Div(id="clustering_parameters", children=[
		    	html.H3(children="Clustering parameters"),
			    
			    html.Div(id='n_neighbors_slider_output_container',
			    		 style={'margin-top': 20}),
			    dcc.Slider(
			        id='n_neighbors_slider',
			        min=1,
			        max=200,
			        step=1,
			        value=20,
			        marks={
			        	1: "1",
			        	5: "5",
			        	10: "10",
			        	20: "20(default)",
			        	50: "50",
			        	125: "125",
			        	200: "200"
		    		},
			    ),
				html.Div([
				    html.P('''
				    	This is roughly related to the minimum number of cells that
						will be grouped into a cluster together, and changes the way
						the UMAP projection is structured.
						''')
					], 
					style={'marginBottom': 20, 'marginTop': 20}
				),

			    html.Div(id='clustering_resolution_slider_output_container',
			    		 style={'margin-top': 20}),			
			    dcc.Slider(
			        id="clustering_resolution_slider",
			        min=0,
			        max=5,
			        step=0.1,
			        value=0.5,
					marks={
			        	0.1: "0.1",
			        	0.25: "0.25",
			        	0.5: "0.5(default)",
			        	1: "1",
			        	2.5: "2.5",
			        	5: "5"
		    		},
			    ),    
				html.Div([
				    html.P('''
				    	Higher clustering resolution leads to a greater number of 
				    	clusters (finer-grained).
						''')
					], 
					style={'marginBottom': 20, 'marginTop': 20}
				),
			],
			style={"display": "none"}),

		    html.Div(children=[
			    html.Button("Load old analysis", 
			    			id="load_analysis_button"),
	    		html.P("Press this button to get started by loading the original analysis and gene lists."),
			    html.Button("Save analysis", 
			    			id="save_analysis_button", 
			    			style={"display": "none"}),
		    ]
		    ),	

		    html.Div(children=[
			    html.Button("Recalculate projection & update plot", 
			    			id="refresh_projection_button"),
			    html.Button("Recalculate clusters & update plot", 
			    			id="refresh_clustering_button"),
			    html.Button("Recalculate everything & update plot", 
			    			 id="refresh_all_button")
		    ],
		    style={'display': 'none'}),

		    html.Div(children=[
			    html.Button("Recalculate pseudotime & update plot", 
			    			id="refresh_pseudotime_button",
			    			style={"display": "none"})
		    ]
		    ),		    

		    html.Div(children=[
			    html.Button("Define new cluster", 
			    			id="define_cluster_button"),
			    html.P("Add currently selected cells from the clustering plot to a new cluster in the user-cluster plot group you have currently selected", id="define_cluster_text")
		    ],
		    style={'marginBottom': 20, 'marginTop': 20}),	
		    

		    # Plots
			html.Div(children=[
				html.Div(children=[
				   	# clustering plot
				    html.H3(children="Clustering plot"),
				    html.P("Show either leiden (automatic) cluster assignments or user-defined cluster assignments"),
				    cc.clustering_dropdown(),
				    cc.plot_clustering_UMAP(),
				    ], className="six columns"
				),
			    html.Div(children=[
				   	# pseudotime plot
				    html.H3(children="Pseudotime plot"),
				    html.P("Show either the pseudotime assignment or calculated differentiation potential of cells, based on the palantir markov-chain-based algorithm"),
				    cc.pseudotime_dropdown(),
				    cc.plot_pseudotime_UMAP(),
				    ], className="six columns"
				)
		    ], className="row"
		    ),

		    html.Div(children=[
				html.Div(children=[
				   	# expression plot
				    html.H3(children="Gene expression (projection)"),
				    html.P("Visualize expression of a single highly-variable gene across single cells, and pull up data from Flybase (Oct, 2019) on the selected gene"),
				    cc.single_gene_dropdown(),
				    cc.plot_expression_UMAP(),
				    html.H3(children="Gene information"),
				    cc.gene_data_table()
				    ], className="six columns"
				),
			    html.Div(children=[
				   	# pseudotime gene expression plot
				    html.H3(children="Gene expression (pseudotime/bulk)"),
				    html.P("Visualize the expression of multiple highly-variable genes against pseudotime or across all cells from the entire dataset"),
				    cc.multi_gene_dropdown(),
				    cc.plot_gene_pseudotime(),
				    cc.plot_gene_violin()
				    ], className="six columns"
				)
		    ], className="row"
		    )

		    ]
		),
	])
