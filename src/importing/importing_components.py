import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

def importing_dataset_dropdown():
    m = dbc.Col(children=[
            dbc.Row(children=[
                dbc.Col(
                    dcc.Dropdown(
                        id='importing_dataset_dropdown',
                        options=[
                            {'label': 'Michki et. al. (2020) fly type-II progenies', 'value': "00001"},
                            {'label': 'Cocanougher et. al. (2020) whole fly CNS', 'value': "00002"},
                            {'label': 'Davie et. al. (2018) aging fly brain', 'value': "00003"},
                            {'label': '10Xv3 5K PBMC', 'value': "00004"}, 
                            {'label': 'Sharma et. al. (2020) mouse somatosensory neurons', 'value': "00005"},
                            {'label': 'Zeisel et. al. (2018) mouse nervous system', 'value': "00006"},
                        ],
                        value=None,
                        placeholder="Select a pre-made dataset",
                        searchable=True,
                        multi=False,
                    ),
                width=9),
                dbc.Col(children=[
                    html.Div(children=[
                        dbc.Button(
                            children=["Load selected dataset"],
                            id="importing_load_dataset_button",
                        ),
                        dcc.Loading(children=[
                            html.Div(children="", 
                                     id='load_selected_dataset_success_output',
                                     style={'margin-top': 10}
                            )
                        ])
                    ])
                ], width=3),
            ]),
        ], width=8)
    return m

def importing_data_upload():
    m = dbc.Col(children=[
            dcc.Upload(
                id='upload_raw_data',
                children=html.Div([
                    '''Drag and drop either an h5ad anndata object,
                    or a .zip file containing your 10X output
                    directory's contents. Alternatively, ''',
                    html.A('click here to select a file.')
                ]),
                style={
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
            html.Div(id='upload_raw_data_success_output',
                     style={'margin-top': 20})
        ], width=4)
    return m