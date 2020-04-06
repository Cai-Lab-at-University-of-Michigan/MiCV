import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from . import processing_components as cc

demo = False 

def processing_layout():
	m = dbc.Tab(label='Processing', children=[
			# File uploading
			dbc.Row(children=[
	    		dbc.Col(children=[
	    			dbc.Card(children=[
	        			dbc.CardHeader(children=[
	        				dbc.Button(
				            "Select data",
				            id="upload-collapse-button",
				            className="mb-3",
				            color="primary",
				        	),
				        ]),
				        dbc.Collapse(children=[
				        	html.Div([
							    html.P('''
							    	Either upload your own data (disabled) 
							    	or search for a pre-made dataset
							    	to load into MiCV for analysis.
									''')
								], 
								style={'marginBottom': 20, 'marginTop': 20}
							),
				            dbc.Row(children=[
				            	cc.processing_dataset_dropdown(),
				            	(cc.processing_data_upload() if (demo is False) else html.Div()),
							], no_gutters=True),
							html.Div(id='upload_spacer_div',
						    		 style={'margin-top': 250}),			
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
							    min=0,
							    max=20000,
							    step=25,
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
							    	unique genes that they must express. Cells expressing very few
							    	genes are often really debris. Cells expressing very many
							    	unique genes are often really doublets (2 cells 1 droplet). However,
							    	depending on the biological makeup of your sample, cells
							    	that express many unique genes might represent
							    	larger-than-average or transcriptionally very active cells.
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
							    	in very few cells might be too poorly recovered to be useful
							    	for downstream analysis, and should be removed. Genes will be further
							    	filtered down even further after this by the
							    	highly variable gene selection process.
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
							    	based on the cell-ranger style of highly-variable gene selection, 
							    	and 2000 is generally a good starting point.
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
							    	The number of neighbors is roughly related to the 
							    	minimum number of cells that
									will be grouped into a cluster together. It changes the way
									the UMAP projection is structured, and can affect cluster
									sizes.
									''')
								], style={'marginBottom': 20, 'marginTop': 20}
							),
							html.Div([
							    html.P('''
							    	The standard neighborhood identification algorithm is 
							    	what you want to use, unless your data has multiple 
							    	biological batches and an obs column labelled "batch".
									''')
								], style={'marginBottom': 10, 'marginTop': 20}
							),
							cc.neighbors_method_radio(),
							dbc.Button("Recalculate projection", 
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
						        step=0.05,
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
							    	clusters (fine-grained), whereas lower clustering resolution
							    	generates fewer cluster (coarse-grained). It is highly encouraged
							    	that you try multiple resolutions. Gauge their fitness by checking
							    	their marker gene overlap (next tab) and through 
							    	in situ validation experiments.
									''')
								], 
								style={'marginBottom': 20, 'marginTop': 20}
							),			    
						    dbc.Button("Recalculate clusters", 
						    			id="refresh_clustering_button"),
						], id="clustering-collapse"),
					]),
	    		],  className="accordion"),
			]),

			dbc.Row(children=[
				dbc.Col(children=[
					# Buttons to do the recalculations
					html.P('''
				    	After loading a dataset and choosing QC,
				    	projection, and clustering parameters, click
				    	this button to prepare your data for analysis.

				    	(Note that this is optional for pre-made
				    	datasets, as basic pre-processing has already
				    	been completed.)
					''', style={'margin-top': 10}),
				    dbc.Button("Recalculate everything", 
				    			 id="refresh_all_button"),
				    html.Div(id='refresh_all_status',
						     style={'margin-top': 10})
		    	], width=6),
			]),

		    dbc.Row(children=[
			    dbc.Col(children=[
				    html.H3(children="Processing UMAP plots"),
				    html.P("Use this dropdown menu to observe your data"),
				    cc.processing_UMAP_dropdown(),
				    cc.n_dims_processing_radio(),
			    	dcc.Loading(children=[cc.plot_processing_UMAP()])	
				], width=6),
				dbc.Col(children=[
				    html.H3(children="Processing QC plots"),
				    html.P("Use this dropdown menu view different QC facets"),
				    cc.processing_QC_dropdown(),
			    	dcc.Loading(children=[cc.plot_processing_QC()])		
				], width=5)
			]),
		])

	return m
