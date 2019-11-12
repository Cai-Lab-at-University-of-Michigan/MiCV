import dash_core_components as dcc
import dash_html_components as html
import custom_components as cc

main_layout = html.Div(children=[
    # title matter
    html.H1(children="MiCV"),
    html.Div(children='''
        A platform for simultaneous cluster-based & pseudotemporal analysis of scRNA-seq data
    '''),

    # clustering parameters
    html.Div(id="analysis_page", children=[
    	html.H3(children="Clustering parameters"),
	    
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
	    html.Div(id='n_neighbors_slider_output_container',
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
	    html.Div(id='clustering_resolution_slider_output_container',
	    		 style={'margin-top': 20}),
	    
	    html.Button("Recalculate & update plot", id="refresh_all_button")
	    ]
	),

    # clustering plot

    html.H3(children="Clustering plot"),
    cc.plot_clustering_UMAP()
])
