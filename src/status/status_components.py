import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

def status_progress():
	m = dbc.Progress(id="status_progress")
	return m

def status_history():
	m = html.Div(id="status_history", children=[""])
	return m

def status_state():
	m = html.Div(id="status_state", children=[""])
	return m
