import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

def exporting_obs_column_dropdown():
    m = dcc.Dropdown(
        id='exporting_obs_column_dropdown',
        options=[
        ],
        value=None,
        multi=False,
        searchable=True,
        placeholder="obs column selection"
        )
    return m

def exporting_obs_value_dropdown():
    m = dcc.Dropdown(
        id='exporting_obs_value_dropdown',
        options=[
        ],
        value=None,
        multi=True,
        searchable=True,
        placeholder="obs column value selection"
        )  
    return m

def exporting_obs_subset_button():
    m = dbc.Button("Generate subset", 
                   id="exporting_obs_button")
    return m

def exporting_obs_subset_link_button():
    m = html.A(dbc.Button("Export subset",
               id="exporting_obs_link_button",
               disabled=True),
               id="exporting_obs_subset_link",
               href=""
        )
    return m

def exporting_button():
    m = html.A(dbc.Button("Export all", 
               id="exporting_button",
               disabled=True),
               id="exporting_link",
               href=""
        )
    return m
