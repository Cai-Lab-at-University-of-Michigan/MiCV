import dash_html_components as html
import dash_bootstrap_components as dbc

def exporting_layout():
	m = dbc.Tab(label="Export/download analysis", children=[
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
			]) # end download analysis tab
	return m