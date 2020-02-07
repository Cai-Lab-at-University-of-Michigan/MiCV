import dash_core_components as dcc
import dash_html_components as html
import custom_components as cc
import uuid

def main_layout(session_id = None):
	if (session_id is None):
		session_id = "testing"
		#session_id = str(uuid.uuid4())

	return session_id, html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    
	    # title matter
	    html.H1(children="MiCV"),
	    html.H3(children="A Multi-Informatic Cellular Visualization tool"),

	    # Links to other parts of MiCV
	    #dcc.Link('Go to processing', href='/apps/processing'),
	    # clustering parameters
	    dcc.Tabs([

	    	### PROCESSING TAB ###
        	
			dcc.Tab(label='Processing', children=[

				# File uploading
				dcc.Upload(
			        id='upload_raw_data',
			        children=html.Div([
			            '''Drag and drop a .zip file containing your 10X output
			            directory's contents, or ''',
			            html.A('select a file')
			        ]),
			        style={
			            'width': '100%',
			            'height': '60px',
			            'lineHeight': '60px',
			            'borderWidth': '1px',
			            'borderStyle': 'dashed',
			            'borderRadius': '5px',
			            'textAlign': 'center',
			            'margin': '10px'
			        },
			        # Allow multiple files to be uploaded
			        multiple=False
			    ),
			    html.Div(id='upload_raw_data_success_output',
				style={'margin-top': 20}),

				# Processing parameters
				html.Div(id="processing_parameters", children=[
			    	html.H3(children="Processing parameters"),
				    
				    html.Div(id='min_max_genes_slider_output_container',
				    		 style={'margin-top': 20}),			
				    dcc.RangeSlider(
				    	id="min_max_genes_slider",
					    min=1,
					    max=20000,
					    step=10,
					    marks={
					        200: "200\n(default min)",
					        1000: "1000",
					        5000: "5000",
					        10000: "10000\n(default max)",
					        15000: "12500"
						},
					    value=[200, 10000]
					),
					html.Div([
					    html.P('''
					    	Filter cells based off of the minimum and maximum number of 
					    	unique genes that they must express. Low/high expressing cells
					    	are often debris or doublets, respectively.
							''')
						], 
						style={'marginBottom': 20, 'marginTop': 20}
					),
					
					html.Div(id='min_cells_slider_output_container',
				    		 style={'margin-top': 20}),
				    dcc.Slider(
				        id="min_cells_slider",
				        min=1,
				        max=200,
				        step=1,
				        value=2,
				        marks={
				        	2: "2\n(default)",
				        	10: "10",
				        	50: "50",
				        	100: "100",
				        	200: "200"
			    		},
				    ),
					html.Div([
					    html.P('''
					    	Filter genes based on how many cells express it. Genes expressed
					    	in fewer than a few cells might be too poorly recovered to be useful
					    	for downstream analysis.
							''')
						], 
						style={'marginBottom': 20, 'marginTop': 20}
					),

					html.Div(id='n_top_genes_slider_output_container',
				    		 style={'margin-top': 20}),
				    dcc.Slider(
				        id="n_top_genes_slider",
				        min=100,
				        max=10000,
				        step=25,
				        value=2000,
				        marks={
				        	500: "500",
				        	1000: "1000",
				        	2000: "2000\n(default)",
				        	5000: "5000",
				        	10000: "10000"
			    		},
				    ),
					html.Div([
					    html.P('''
					    	Choose how many highly variable genes to consider for downstream
					    	dimensionality reduction, clustering, and other analysis. This is
					    	based on the cell-ranger style of highly-variable gene selection.
							''')
						], 
						style={'marginBottom': 20, 'marginTop': 20}
					),

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
				        	20: "20\n(default)",
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
				        	0.5: "0.5\n(default)",
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
				]),

				# Buttons to do the recalculations
			    html.Div(children=[
				    html.Button("Recalculate projection", 
				    			id="refresh_projection_button"),
				    html.Button("Recalculate clusters", 
				    			id="refresh_clustering_button"),
				    html.Button("Recalculate everything", 
				    			 id="refresh_all_button")
			    ]),

			    html.Div(children=[
				    html.Button("Recalculate pseudotime", 
				    			id="refresh_pseudotime_button")
			    ]),
			    html.Div(children=[
				    html.Div(children=[
					   	# pseudotime gene expression plot
					    html.H3(children="Processing UMAP plots"),
					    html.P("Use this dropdown menu to observe your data"),
					    cc.processing_UMAP_dropdown(),
				    	cc.plot_processing_UMAP()		
					    ], className="six columns"),
					html.Div(children=[
						   	# pseudotime gene expression plot
						    html.H3(children="Processing QC plots"),
						    html.P("Use this dropdown menu view different QC facets"),
						    cc.processing_QC_dropdown(),
					    	cc.plot_processing_QC()		
						    ], className="six columns"
					),
					], className="row"
				)
		    ]),

		    ### ANNOTATION TAB ###
			dcc.Tab(label='Annotation', children=[

			    html.Div(children=[
				    html.Button("Load old analysis", 
				    			id="load_analysis_button"),
		    		html.P("Press this button to get started by loading the original analysis and gene lists."),
				    html.Button("Save analysis", 
				    			id="save_analysis_button", 
				    			style={"display": "none"}),
				]),  

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
	])