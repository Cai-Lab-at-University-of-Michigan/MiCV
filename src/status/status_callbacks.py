import dash
from dash.dependencies import Input, Output, State

from helper_functions import *

from app import app

from . status_functions import *


@app.callback(
    [Output("status_progress", "value"), 
     Output("status_progress", "children")],
    [Input("status_interval", "n_intervals")],
    [State("session-id", "children"),
     State("status_progress", "value")]
)
def update_progress(n, session_ID, previous_progress):
    default_return = [dash.no_update, dash.no_update]

    progress = cache_progress(session_ID)
    if ((n in ["", None, 0, []])
    or	(progress in ["", None, []])
    or	(progress == previous_progress)):
    	return default_return
    else:
    	return [progress, str(progress) + " %"]


@app.callback(
    Output("status_history", "children"),
    [Input("status_interval", "n_intervals")],
    [State("session-id", "children"),
     State("status_history", "children")]
)
def update_history(n, session_ID, previous_history):
    default_return = dash.no_update

    history = build_history_table(session_ID)
    if ((n in ["", None, 0, []])
    or	(history in ["", None, []])):
    	return default_return
    else:
    	return history


@app.callback(
    Output("status_state", "children"),
    [Input("status_interval", "n_intervals")],
    [State("session-id", "children"),
     State("status_state", "children")]
)
def update_state(n, session_ID, previous_state):
    default_return = dash.no_update

    if (n in ["", None, 0, []]):
    	return default_return
    
    state = build_state_table(session_ID)
    if (state in ["", None, []]):
    	return default_return
    else:
    	return state