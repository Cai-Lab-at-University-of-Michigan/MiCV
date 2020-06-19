import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

from . import exporting_components as cc

def exporting_layout(demo=False):
	if (demo is False):
		m = dbc.Tab(label="Export/download analysis", children=[
					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Download analysis"),
							html.P(
								'''
								Download a copy of your entire scRNA-seq data.
								This data is a scanpy-readable anndata object, ready
								for further analysis using scanpy.
								'''),
							cc.exporting_button()
						], width=6),
					]),

					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Export subset"),
							html.P(
								'''
								Download a copy of a SUBSET of your scRNA-seq data. Choose
								an observation column to subset by, then choose which categories
								from that column you'd like to keep.
								This file will be ready for re-uploading to MiCV or scanpy for
								your next iteration of analysis.
								'''),
							dbc.Row(children=[
								dbc.Col(children=[
									cc.exporting_obs_column_dropdown()
								], width=4),
								dbc.Col(children=[
									cc.exporting_obs_value_dropdown()
								], width=4),
								dbc.Col(children=[
									cc.exporting_obs_subset_button()
								], width=2),
								dbc.Col(children=[
									cc.exporting_obs_subset_link_button()
								], width=2)
							]),
						], width=6),
					]),
				]) # end download analysis tab
	else:
		m = dbc.Tab(label="Export/download analysis", children=[
					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Download analysis"),
							html.P(
								'''
								Download functionality is disabled,
								as this is a demo instance of MiCV.
								'''),
						], width=6),
					])
				])
	return m