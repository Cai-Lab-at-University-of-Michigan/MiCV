import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

from . importing_functions import *

def importing_dataset_dropdown():
    m = dbc.Card(children=[
            
            dbc.CardHeader(children=[
                html.H4("Sample datasets", id="sample_datasets_title"),
                dbc.Tooltip(
                    '''
                    Use this interface to load sample, openly-accessible datasets. 
                    Some have been downsampled (fewer cells) from their respective
                    publication's reported size.
                    ''',
                    target="sample_datasets_title",
                    placement="above",
                    delay={"show": 50, "hide": 100}
                ),
            ]),
            
            dbc.CardBody(children=[
                dbc.Container(children=[
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
                        width=10),
                        dbc.Col(children=[
                            dbc.Button(
                                children=["Load"],
                                id="importing_load_dataset_button",
                            ),
                        ], width=2),
                    ])
                ], fluid=True)
            ]),

            dbc.CardFooter(children=[
                dcc.Loading(children=[
                    html.Div(children="", 
                             id='load_selected_dataset_success_output',
                             style={'margin-top': 10}
                    )
                ])
            ])
        ])
    return m

def importing_data_upload(demo=False):
    # if it's a demo instance, don't allow uploads
    m = dbc.Card(children=[
        
        dbc.CardHeader(children=[
            html.H4("Upload"),
        ]),

        dbc.CardBody(children=[
            html.P(
            '''
            Data upload disabled
            '''
            )
        ])
    ], style={"background-color": "rgba(0, 0, 0, 0.15)"})
    
    if (demo is True):    
        return m

    from flask_security import current_user
    try:
        email = current_user.email
    except:
        return m

    m = dbc.Card(children=[
            
            dbc.CardHeader(children=[
                html.H4("Upload", id="upload_title"),
                dbc.Tooltip(
                    '''
                    Currently, only the following datatypes are accepted
                    for importing into MiCV:

                    an h5ad file from MiCV or scanpy, or

                    a zipped directory, containing the
                    output of 10X Cell Ranger (-style) mapping
                    ''',
                    target="upload_title",
                    placement="above",
                    delay={"show": 50, "hide": 100}
                ),
            ]),
            
            dbc.CardBody(children=[
                dbc.Container(children=[
                    dcc.Upload(
                        id='upload_raw_data',
                        children=html.Div([
                            html.A(
                                '''
                                Drag and drop your file, 
                                or click here to select a 
                                file for upload.
                                '''
                            )
                        ], id="upload_raw_data_div"),
                        style={
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '20px'
                        },
                        multiple=False
                    ),
                ])
            ]),

            dbc.CardFooter(children=[
                dcc.Loading(children=[
                    html.Div(children="", 
                             id='upload_raw_data_success_output',
                             style={'margin-top': 10}
                    )
                ])
            ])
        ])
    return m

def importing_user_dataset_list(demo=False):
    # if it's a demo instance, there will be no files to show
    m = dbc.Card(children=[
        
        dbc.CardHeader(children=[
            html.H4("Your saved datasets", id="saved_datasets_title")
        ]),

        dbc.CardBody(children=[
            html.P(
            '''
            Dataset saving disabled
            '''
            )
        ])
    ], style={"background-color": "rgba(0, 0, 0, 0.15)"})
    
    if (demo is True):
        return m
    
    from flask_security import current_user
    try:
        email = current_user.email
    except:
        return m
    try:
        dataset_list = get_dataset_list(email)
    except:
        dataset_list = []
    if (dataset_list in [[], None, ""]):
        c = dbc.RadioItems(
            options=[{"label": "You have no datasets saved", 
                      "value": "no_datasets_saved"}],
            value="no_datasets_saved",
            id="importing_user_dataset_list_radio",
        )
        disabled = True
    else:
        c = dbc.RadioItems(
                options=[{"label": d[0], "value": d[0]} for d in dataset_list],
                value=dataset_list[0][0],
                id="importing_user_dataset_list_radio",
            )
        disabled = False

    d = dbc.Row(children=[
            dbc.Col(children=[c], width=9),
            dbc.Col(children=[
                dbc.ButtonGroup(children=[
                    dbc.Button("Load", 
                               id="importing_user_dataset_load_button",
                               color="primary",
                               disabled=disabled),
                    #dbc.Button("Share", 
                    #           id="importing_user_dataset_share_button",
                    #           color="secondary"),
                    dbc.Button("Delete", 
                               id="importing_user_dataset_delete_button",
                               color="danger",
                               disabled=disabled),
                ], vertical=True)
            ], width=3)
        ])
    m = dbc.Card(children=[
            
        dbc.CardHeader(children=[
            html.H4("Your saved datasets", id="saved_datasets_title"),
            dbc.Tooltip(
                '''
                Use this interface to load or delete previously saved
                datasets. This will only show datasets saved internally
                through MiCV's save/export tab, not datasets that were
                exported to your computer. Your most recently saved 
                datasets are listed at the top.
                ''',
                target="saved_datasets_title",
                placement="above",
                delay={"show": 200, "hide": 200}
            ),
        ]),
            
            dbc.CardBody(children=[
                d
            ]),
            
            dbc.CardFooter(children=[
                dcc.Loading(children=[
                    html.Div(children="", 
                             id='importing_user_dataset_success_output',
                             style={'margin-top': 10}
                    )
                ])
            ])
        ])
    return m

def importing_greeting(demo=False):
    m = dbc.Jumbotron(children=[
            html.H3("Welcome to the MiCV demo"),
            html.P("Let's get started.")
        ])

    if (demo is False):
        from flask_security import current_user
        try:
            email = current_user.email
        except:
            return m
        
        m = dbc.Jumbotron(children=[
            html.H3("Welcome back, " 
                    + str(current_user.email).split("@",1)[0]),
            html.P("Let's get started.")
        ])

    return m
    
