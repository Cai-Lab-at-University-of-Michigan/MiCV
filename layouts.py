import dash_core_components as dcc
import dash_html_components as html
import custom_components as cc
import uuid

def main_layout():
	session_id = str(uuid.uuid4())
	return html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    # title matter
	    html.H1(children="MiCV"),
	    html.Div(children='''
	        A platform for simultaneous cluster-based & pseudotemporal analysis of scRNA-seq data
	    '''),

	    # clustering parameters
	    html.Div(id="analysis_page", children=[
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
		    
		    html.Div(children=[
			    html.Button("Recalculate projection & update plot", 
			    			id="refresh_projection_button"),
			    html.Button("Recalculate clusters & update plot", 
			    			id="refresh_clustering_button"),
			    html.Button("Recalculate everything & update plot", 
			    			 id="refresh_all_button")
		    ]
		    ),
		    

		    # Plots
			html.Div(children=[
				html.Div(children=[
				   	# clustering plot
				    html.H3(children="Clustering plot"),
				    cc.plot_clustering_UMAP(),
				    ], className="six columns"
				),
			    html.Div(children=[
				   	# pseudotime plot
				    html.H3(children="Pseudotime plot"),
				    cc.plot_pseudotime_UMAP(),
				    ], className="six columns"
				)
		    ], className="row"
		    )

		    ]
		),
	])
