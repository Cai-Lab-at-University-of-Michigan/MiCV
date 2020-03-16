import dash
from dash.dependencies import Input, Output, State

from helper_functions import *
from app import app

@app.callback(
	Output("null_container_2", "children"),
	[Input("cleanup_interval", "n_intervals")]
)
def run_cleanup(n_intervals):
	if not (n_intervals in [None, 0, "", []]):
		remove_old_cache(n_days=1.5)

	return dash.no_update