import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

from . import exporting_components as cc

def exporting_layout(demo=False):
	m = dbc.Tab(label="Save/export analysis", children=[
					dbc.Row(children=[
						dbc.Col(children=[
							html.H3(children="Export analysis"),
							html.P(
								'''
								Download functionality is disabled.
								'''),
						], width=6, 
						style={"background-color": "rgba(0, 0, 0, 0.15)"}),
					])
				])

	if (demo is True):
		return m
	from flask_security import current_user
	try:
		email = current_user.email
	except:
		return m

	m = dbc.Tab(label="Save/export analysis", children=[
				dbc.Row(children=[
					dbc.Col(children=[
						html.H3(children="Save analysis"),
						html.P(
							'''
							Save this analyzed dataset to your MiCV 
							account, to come back to later.
							'''),
						dbc.Row(children=[
							dbc.Col(children=[
								cc.save_dataset_input()], width=4),
							dbc.Col(children=[
								cc.save_dataset_button()], width=4)
						])
					], width=8)
				]),

				dbc.Row(children=[
					dbc.Col(children=[
						html.H3(children="Export analysis"),
						html.P(
							'''
							Download a copy of your scRNA-seq data. 
							Choose either to take all cells (default)
							or a subset of cells that you choose using
							the dropdown menus here.
							This file will be ready for re-uploading 
							to MiCV or scanpy for further analysis.
							'''),
						dbc.Row(children=[
							dbc.Col(children=[
								cc.exporting_subset_radio()
							], width=4)
						]),
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
		
	return m