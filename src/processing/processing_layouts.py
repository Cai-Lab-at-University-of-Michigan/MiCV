import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from . import processing_components as cc

def processing_layout():
	m = dbc.Tab(label='Processing', children=[
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
							html.Div(id='n_dims_proj_radio_output_container',
						    		 style={'margin-top': 20}),
							cc.n_dims_proj_radio(),
							html.Div([
							    html.P('''
							    	How many dimensions for the UMAP projection - choose
							    	3 if you'd like to have the ability to do 
							    	3D projections, at the expense of having slightly 
							    	compressed 2D projections.
									''')
								], style={'marginBottom': 20, 'marginTop': 20}
							),
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
				    html.Div(id='refresh_all_status',
						     style={'margin-top': 20})
			    ], width=3)
			]),

		    dbc.Row(children=[
			    dbc.Col(children=[
				   	# pseudotime gene expression plot
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
				], width=6)
			]),
		])

	return m