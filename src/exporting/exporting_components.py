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
        placeholder="(obs) column to subset by"
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
        placeholder="(obs) column values to keep"
        )  
    return m

def exporting_obs_subset_button():
    m = dbc.Button("Generate file", 
                   id="exporting_obs_button")
    return m

def exporting_obs_subset_link_button():
    m = html.A(dbc.Button("Download file",
               id="exporting_obs_link_button",
               disabled=True),
               id="exporting_obs_subset_link",
               href=""
        )
    return m
'''
def exporting_button():
    m = html.A(dbc.Button("Export all", 
               id="exporting_button",
               disabled=True),
               id="exporting_link",
               href=""
        )
    return m
'''
def exporting_subset_radio():
    m = dbc.RadioItems(options=[
                {"label": "All cells (default)", "value": "all"},
                {"label": "Subset of cells", "value": "subset"}
            ], 
            id="exporting_subset_radio",
            value="all"
        )
    return m

def save_dataset_input():
    m = dbc.FormGroup(children=[
            dbc.Input(id="save_dataset_name", 
                      placeholder="Enter a save name for this dataset",
                      type="text",
                      valid=None),
            dbc.FormFeedback("Data saved", 
                              id="save_dataset_success_feedback",
                              valid=True),
            dbc.FormFeedback("Error: data not saved", 
                              id="save_dataset_failure_feedback",
                              valid=False),
          ], id="save_dataset_name_formgroup")
    return m

def save_dataset_button():
    m = dbc.Row(children=[
            dbc.Col(children=[
                dbc.Button("Save dataset to MiCV",
                            id="save_dataset_button"),
                dcc.Loading(children=[
                    html.Div(children=[" "],
                             id="save_dataset_output"
                    )
                ])
            ])
        ])
    return m
