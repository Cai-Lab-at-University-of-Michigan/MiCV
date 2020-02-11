import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import custom_components as cc
import uuid

def main_layout(session_id = None):
	if (session_id is None):
		session_id = "testing"
		#session_id = str(uuid.uuid4())

	return session_id, html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_1', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI

	    dbc.Container(fluid=True, children=[
		    # title matter
		    html.H1(children="MiCV"),
		    html.H3(children="A Multi-Informatic Cellular Visualization tool"),

		    # Links to other parts of MiCV
		    #dcc.Link('Go to processing', href='/apps/processing'),
		    # clustering parameters
		    dbc.Tabs([

		    	### PROCESSING TAB ###
	        	
				dbc.Tab(label='Processing', children=[

					# File uploading
					dbc.Col(children=[
						html.Div(children=[
		        			dbc.Card(children=[
			        			dbc.CardHeader(children=[
			        				dbc.Button(
						            "Upload",
						            id="upload-collapse-button",
						            className="mb-3",
						            color="primary",
						        	),
						        ]),
						        dbc.Collapse(children=[
						            dcc.Upload(
								        id='upload_raw_data',
								        children=html.Div([
								            '''Drag and drop either an h5ad anndata object,
								            or a .zip file containing your 10X output
								            directory's contents. Alternatively ''',
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
								], id="upload-collapse"),
							]),

		        			dbc.Card(children=[
			        			dbc.CardHeader(children=[
			        				dbc.Button(
						            "QC parameters",
						            id="QC-collapse-button",
						            className="mb-3",
						            color="primary",
						        	),
						        ]),
						        dbc.Collapse(children=[
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
										], style={'marginBottom': 20, 'marginTop': 20}
									),
								], id="QC-collapse"),
						    ]),

							dbc.Card(children=[
			        			dbc.CardHeader(children=[
			        				dbc.Button(
						            "Projection parameters",
						            id="projection-collapse-button",
						            className="mb-3",
						            color="primary",
						        	),
						        ]),

						        dbc.Collapse(children=[
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
										], style={'marginBottom': 20, 'marginTop': 20}
									),
									cc.neighbors_method_radio(),
									dbc.Button("Recalculate (only) projection", 
							    				id="refresh_projection_button"),
								], id="projection-collapse"),
						    ]),

		        			dbc.Card(children=[
			        			dbc.CardHeader(children=[
			        				dbc.Button(
						            "Clustering parameters",
						            id="clustering-collapse-button",
						            className="mb-3",
						            color="primary",
						        	),
						        ]),
						        dbc.Collapse(children=[
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
								    dbc.Button("Recalculate clusters", 
								    			id="refresh_clustering_button"),
								], id="clustering-collapse"),
						    ]),
		        		],  className="accordion"),
						dbc.Col(children=[
							# Buttons to do the recalculations
						    dbc.Button("Recalculate everything", 
						    			 id="refresh_all_button"),
						    html.P('''
						    	Run this first after uploading your data 
						    	& doing QC parameter selection
								'''),
						    dbc.Button("Recalculate pseudotime", 
						    			id="refresh_pseudotime_button"),
						    html.P('''
						    	Be sure to select a starter cell in the plot 
						    	below before running this.\nThis will take a
						    	lonnnngggg time and might fail to converge; please
						    	be patient.
								'''),
					    ], width=3)
					]),
			
				    dbc.Row(children=[
					    dbc.Col(children=[
						   	# pseudotime gene expression plot
						    html.H3(children="Processing UMAP plots"),
						    html.P("Use this dropdown menu to observe your data"),
						    cc.processing_UMAP_dropdown(),
					    	dcc.Loading(children=[cc.plot_processing_UMAP()])		
						], width=6),
						dbc.Col(children=[
							   	# pseudotime gene expression plot
							    html.H3(children="Processing QC plots"),
							    html.P("Use this dropdown menu view different QC facets"),
							    cc.processing_QC_dropdown(),
						    	cc.plot_processing_QC()		
						], width=6)
					]),
				]),


				### MARKER GENE TAB ###
				dbc.Tab(label="Marker genes", children=[
					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Clustering plot"),
							html.P(
								'''
								obs columns with fewer than 50 unique values
								can be used to define groups for marker gene 
								identification (typically you should use the 
								leiden column or one of the user_n manually 
								annotated columns)
								'''),
							cc.marker_gene_UMAP_dropdown(),
							cc.plot_marker_gene_UMAP()
						], width=6),
					]),
					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Marker gene identification"),
							html.P(
								'''
								Select which groups you would like to identify
								marker genes for and which method you would like
								to use for identification
								'''),
							cc.marker_gene_group_dropdown(),
							cc.marker_gene_method_radio(),
							dbc.Button("Recalculate marker genes", 
						    			id="recalculate_marker_genes"),
							cc.marker_gene_plot()
						]),
					])

				]), # end marker gene tab
				
			    ### ANNOTATION TAB ###
				dbc.Tab(label='Annotation', children=[
				    html.Div(children=[
					    dbc.Button("Load processed data", 
					    			id="load_analysis_button"),
			    		html.P("Press this button to get started by loading the processed data"),
					    html.Button("Save analysis", 
					    			id="save_analysis_button", 
					    			style={"display": "none"}),
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
						    ], width=6,
						),
					    dbc.Col(children=[
						   	# pseudotime plot
						    html.H3(children="Pseudotime plot"),
						    html.P("Show either the pseudotime assignment or calculated differentiation potential of cells, based on the palantir markov-chain-based algorithm"),
						    cc.pseudotime_dropdown(),
						    cc.plot_pseudotime_UMAP(),
						    ], width=6,
						)
				    ], 
				    ),

				    dbc.Row(children=[
						dbc.Col(children=[
						   	# expression plot
						    html.H3(children="Gene expression (projection)"),
						    html.P("Visualize expression of a single highly-variable gene across single cells, and pull up data from Flybase (Oct, 2019) on the selected gene"),
						    cc.single_gene_dropdown(),
						    cc.plot_expression_UMAP(),
						    html.H3(children="Gene information"),
						    cc.gene_data_table()
						    ], width=6,
						),
					    dbc.Col(children=[
						   	# pseudotime gene expression plot
						    html.H3(children="Gene expression (pseudotime/bulk)"),
						    html.P("Visualize the expression of multiple highly-variable genes against pseudotime or across all cells from the entire dataset"),
						    cc.multi_gene_dropdown(),
						    cc.plot_gene_pseudotime(),
						    cc.plot_gene_violin()
						    ], width=6,
						)
				    ], 
				    )
				]), ## end annotation tab

				### DOWNLOAD ANALYSIS TAB ###
				dbc.Tab(label="Download analysis", children=[
					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Download analysis"),
							html.P(
								'''
								Download a copy of your scRNA-seq data in h5ad format.
								This data is a scanpy-readable anndata object, ready
								for further analysis using scanpy.
								'''),
							html.A(id='download_anndata_h5ad_button', 
								   children="Download h5ad file",
								   href="/MiCV/download/h5ad") # TODO: remove hardlink
						], width=6),
					]),

				]), # end download analysis tab
			])
		])
])