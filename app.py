from dash import Dash, html, Input, Output, State, callback, dcc, dash_table, callback_context
import dash
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import sqlite3
import pandas as pd
import numpy as np
import os 


###############################################################################################
# Functions
###############################################################################################

##### load data

db_path = os.path.join(os.path.dirname(__file__), "Predictables.db")

def load_master_data():
    conn = sqlite3.connect(db_path)
    df = pd.read_sql("SELECT * FROM Predictables", conn)
    conn.close()
    return df

def load_sim():
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM SIM", conn)
    conn.close()
    return df

##### Search for RI and fragment masses
def get_ri_list(substances, phase):
    col = "RI polar" if phase == "RI polar" else "RI semi-standard non-polar"
    return [s[col] for s in substances if s[col]]

def get_sim_for_cas(cas_list, highest_intensity_only=False):
    df_S = load_sim()
    filtered = df_S[df_S["CAS-No."].isin(cas_list)]
    if highest_intensity_only:
        filtered = filtered[filtered["rel-Intensity"] == 1]
    sim_dict = filtered.groupby("CAS-No.")["Fragment"].apply(lambda x: ", ".join(sorted(map(str, x)))).to_dict()
    return sim_dict

##### Chromatogram
def show_chrom(ri_list, sigma=7):
    fig = go.Figure()

    # initial empty chromatogram
    fig.add_trace(go.Scatter(
        x=[500, 3500],
        y=[0, 0],
        mode='lines',
        line=dict(color='lightgray', width=2)
    ))
    fig.update_layout(
        xaxis=dict(range=[500, 3500], title='Retention index (RI)'),
        yaxis=dict(range=[0, 1.2], title='Peak or no Peak'),
        title='Preview chromatogramm',
        showlegend=False
    )

    # add a peak:
    for i in ri_list:

        # peaks bell-shaped
        x_peak = np.linspace(i-4*sigma, i+4*sigma, 200)
        y_peak = 1 * np.exp(-((x_peak - i) ** 2) / (2 * sigma ** 2))

        fig.add_trace(go.Scatter(
            x=x_peak,
            y=y_peak,
            mode='lines',
            line=dict(width=2)
        ))

        fig.add_annotation(
            x=i,
            y=1,
            text=str(int(i)),
            showarrow=False,
            textangle=270,
            font=dict(size=14, color='black'),
            xanchor='center',
            yanchor='bottom'
        )

    return fig

##### Mass spectrum
def show_ms(cas, title):
    df_sim = load_sim()
    df_sub = df_sim[df_sim["CAS-No."] == cas].sort_values("Fragment")
    
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df_sub["Fragment"],
        y=df_sub["rel-Intensity"],
        marker_color='darkblue',
        width=0.8,
        text=df_sub["Fragment"],
        textposition="outside"
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="m/z",
        yaxis_title="Relative Intensity",
        template="plotly_white"
    )
    return fig

###############################################################################################
# Application
###############################################################################################

app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME]
)  
server = app.server

df = load_master_data()

app.layout = dbc.Container(
    [  
    
    # Overhead Navigation Bar with the Logo and a Title
        dbc.Navbar(
            [
                dbc.NavbarBrand(
                    [
                        html.Img(
                            src=dash.get_asset_url("TUD_Logo.png"),
                            height="80px",
                            style={
                                "marginRight": "20px",
                                "marginLeft": "10px",
                            },
                        ),
                        "PaperOdourPredictables - A tool to improve your GC-MS/O analysis for paper and cardboard",
                    ],
                    href="/",
                    className="mr-4",
                    style={"color": "white"},
                ),
            ],
            color = "#00305d",
        ),

    # introduction  
        dbc.Container(
            [
                html.Div(
                    [
                        html.Br(),
                        html.P([
                            "PaperOdourPredictables is a tool which contains data specific to substances typically found in paper. ",
                            "It can be used for targeting and improving your gas chromatography - mass spectrometry with olfactometry (GC-MS/O) analysis and to interpret and assess the likely sources of odour-active substances in the paper production process. "
                            "The database was built from a broad literature review combined with data from the working group. ",
                            html.Br(),
                            html.Br(),
                            "The odorants typically found in paper and board are listed with data on gas chromatography (GC), mass spectrometry (MS), olfactometry (O), and the likely source in the production process:",
                        ]),
                        html.Ul(
                            [
                                html.Li([
                                    html.Strong("Odour-describing attributes"),
                                    html.Br(),
                                    "from our in-house odour database and from the reference literature"
                                ]),
                                html.Li([
                                    html.Strong("Odour thresholds"),
                                    " for many substances",
                                    html.Br(),
                                    "calculated as the median of the odour detection thresholds in Air (if it is not provided, then another medium) reported by van Gemert, L. J. (2011). Odour Thresholds: Compilations of odour threshold values in air, water and other media (2nd enlarged and revised edition). Oliemans Punter & Partners BV."
                                ]),
                                html.Li([
                                    html.Strong("Retention indices (RI)"),
                                    "  for a polar column (e.g. WAX) and a semi-standard non-polar column (e.g. DB-5)",
                                    html.Br(),
                                    "from NIST Mass Spectrometry Data Center, William E. Wallace, director. 'Retention Indices' in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J. Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, ",
                                    html.A("https://doi.org/10.18434/T4D303", href="https://doi.org/10.18434/T4D303", target="_blank"),
                                    ", (retrieved August 2026)"                                
                                ]),
                                html.Li([
                                    html.Strong("Mass spectra information"),
                                    " (typical fragments and relative intensities)",
                                    html.Br(),
                                    " from NIST Mass Spectrometry Data Center, William E. Wallace, director. 'Mass Spectra' in NIST Chemistry WebBook, NIST Standard Reference Database Number 69, Eds. P.J. Linstrom and W.G. Mallard, National Institute of Standards and Technology, Gaithersburg MD, 20899, ",
                                    html.A("https://doi.org/10.18434/T4D303", href="https://doi.org/10.18434/T4D303", target="_blank"),
                                    ", (retrieved August 2026)" 
                                ]),
                                html.Li([
                                    html.Strong("Origin / Source "),
                                    "in the production process of paper and cardboard ",
                                    html.Br(),
                                    "This list is not exhaustive. In addition to the information listed, other sources are not excluded. ",
                                    html.Br(),
                                    "For each substance, the reference is given in which the substance was described as an odorant derived from paper, as well as, (where applicable) a reference indicating the substance's source or origin. ",
                                    "The bibliography can be downloaded here: ",
                                    html.Br(),
                                    html.Button(
                                        "Download Bibliography.txt",
                                        id="btn-download-bib",
                                        n_clicks=0, 
                                        className="btn btn-primary",
                                        style={
                                            "background-color": "#005d39",
                                            'opacity': '1', 
                                            'border-color': 'transparent',
                                            "width": "250px",
                                            "marginLeft": "1px"
                                        },
                                    ),
                                    dcc.Download(id="download-bib"),
                                ])


                            ],
                            style={
                                'margin-left': '10px',
                                'margin-bottom': '4px'
                            }
                        ),
                        html.Br(),
                        html.Strong("How do I use this application?"),
                        html.Ol(
                            [
                                html.Li("Use the checkboxes to select substances you want to analyse as targets. Then click the button to add the selected substances to the output table."),
                                html.Ul(
                                    [
                                        html.Li("You can also filter the columns using the filter fields below the header. To enable upper and lower case, press the red Aa field on the right-hand side of the filter cells. Enter the text in the filter cell and press Enter."),
                                        html.Li("The number of entries in the table is also displayed."),
                                    ],
                                    style={
                                        'margin-left': '10px',
                                        'margin-bottom': '4px'
                                    }
                                ),
                                html.Li("You will find your selected substances in the output table."),
                                html.Ul(
                                    [
                                        html.Li("To remove substances from this table, select them using the checkboxes and press the remove button."),
                                        html.Li("The columns in this table can also be filtered."),
                                        html.Li("You can download the output table in Excel (.xlsx) or CSV (.csv) format.")
                                    ],
                                    style={
                                        'margin-left': '10px',
                                        'margin-bottom': '4px'
                                    }
                                ),
                                html.Li("The entries in the output table can then be visualised."),
                                html.Ul(
                                    [
                                        html.Li("First, select the phase of your GC column."),
                                        html.Li("Then click the display button."),
                                        html.Li("Now you see a chromatogram with peaks of the selected substances at the corresponding retention indices."),
                                        html.Li("You can also click on the peaks to display the mass spectrum with typical fragments."),
                                    ],
                                    style={
                                        'margin-left': '10px',
                                        'margin-bottom': '4px'
                                    }
                                ),
                            ]
                        ),
                        html.Hr(style={"border-top": "5px solid #00305d"}),
                    ]
                )
            ]
        ),

    # tables
        dbc.Container(
            [
                html.Div(
                    [
                        html.Strong("1) Master table"),
                        html.P("The table contains odour-active substances typically found in paper and board."),
                        html.Div(id="table-count"),
                        dcc.Store(id="selected-substances"),
                        dash_table.DataTable(
                            id="table-master",
                            columns=[{"name": c, "id": c, "selectable": True} for c in ["Odorant", "CAS-No.", "Odour description", "Threshold mg/m³ Air", "RI polar", "RI semi-standard non-polar", "Origin/Source", "Literature paper samples", "Literature origin"]],
                            data=df[["Odorant", "CAS-No.", "Odour description", "Threshold mg/m³ Air","RI polar", "RI semi-standard non-polar", "Origin/Source", "Literature paper samples", "Literature origin"]].to_dict("records"),
                            filter_action="native",
                            row_selectable="multi",
                            selected_rows=[],
                            page_size=10,
                            style_cell={
                                'padding': '5px',
                                'textAlign': 'left',
                                'minWidth': '30px',
                                'width': '150px',
                                'maxWidth': '400px',
                                'whiteSpace': 'normal'
                            },
                            style_header={
                                'backgroundColor': "#00305d",
                                'fontSize' : '18px',
                                'fontWeight': 'bold',
                                'color' : 'white'
                            },
                            style_filter={
                                'backgroundColor': "#dbf1ffd0",
                                'color': 'black',
                            }               
                        ),
                        html.Br(),
                        html.Button(
                            "Add to table",
                            id="btn-add",
                            n_clicks=0, 
                            className="btn btn-primary",
                            style={
                                "background-color": "#00305d",
                                'opacity': '1', 
                                'border-color': 'transparent',
                                "width": "250px",
                            },
                        ),
                        html.Hr(style={"border-top": "5px solid #00305d"}),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("2) Output table"),
                        html.P("Here you can see the selected substances:"),
                        dash_table.DataTable(
                            id="table-selected",
                            columns=[{"name": c, "id": c} for c in ["Odorant", "CAS-No.", "Odour description", "RI polar", "RI semi-standard non-polar", "Fragments"]],
                            data=[],
                            filter_action="native",
                            row_selectable="multi",
                            page_size=10,
                            style_cell={
                                'padding': '5px',
                                'textAlign': 'left',
                                'minWidth': '30px',
                                'width': '150px',
                                'maxWidth': '400px',
                                'whiteSpace': 'normal'
                            },
                            style_header={
                                'backgroundColor': "#243F91",
                                'fontSize' : '18px',
                                'fontWeight': 'bold',
                                'color' : 'white'
                            },
                            style_filter={
                                'backgroundColor': "#dbf1ffd0",
                                'color': 'black',
                            } 
                        ),
                        html.Br(),
                        html.Div(id="output-sim"),
                        html.Div(id="output-sim-1"),
                        html.Br(),
                        html.Button(
                            "Remove from table",
                            id="btn-remove",
                            n_clicks=0, 
                            className="btn btn-primary",
                            style={
                                "background-color": "#00305d",
                                'opacity': '1', 
                                'border-color': 'transparent',
                                "width": "250px",
                            },
                        ),
                        html.Hr(),
                        html.Button(
                            "Download Excel",
                            id="btn-download-excel",
                            n_clicks=0, 
                            className="btn btn-primary",
                            style={
                                "background-color": "#005d39",
                                'opacity': '1', 
                                'border-color': 'transparent',
                                "width": "250px",
                            },
                        ),
                        dcc.Download(id="download-excel"),

                        html.Button(
                            "Download CSV",
                            id="btn-download-csv",
                            n_clicks=0, 
                            className="btn btn-primary",
                            style={
                                "background-color": "#005d39",
                                'opacity': '1', 
                                'border-color': 'transparent',
                                "width": "250px",
                                "marginLeft": "10px"
                            },
                        ),
                        dcc.Download(id="download-csv"),

                        html.Hr(style={"border-top": "5px solid #00305d"}),
                    ]
                )
            ]
        ),
    
    # plots
        dbc.Container(
            [
                html.Div(
                    [
                        html.Strong("3) Visualization"),
                        html.P("The chromatogram and mass spectra are displayed here. First, the phase of the GC column must be selected:"),
                        dcc.Dropdown(
                            id="dropdown-phase",
                            options=[
                                {"label": "polar", "value": "RI polar"},
                                {"label": "semi-standard non-polar", "value": "RI semi-standard non-polar"}
                            ],
                            placeholder="Select the phase",
                            style={"width": "250px"}
                        ),
                        html.Div(id="phase-warning", style={"color": "red"}),
                        html.Br(),
                        html.Button(
                            "Display",
                            id="btn-analyse",
                            n_clicks=0, 
                            className="btn btn-primary",
                            style={
                                "background-color": "#00305d",
                                'opacity': '1', 
                                'border-color': 'transparent',
                                "width": "250px",
                            },
                        ),
                        html.Br(),
                        html.Hr(),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("Chromatogram"),
                        html.P(
                            [
                                "This is an interactive image. Use the left mouse button to select an area and enlarge it by drawing a window. Double-click to reset zoom.",
                                html.Br(),
                                "Click on a peak to see the mass spectrum of the substance."
                            ]
                        ),
                        html.Div(
                            dcc.Graph(
                                id="chromatogram", 
                                figure=show_chrom([]),
                            )
                        ),
                        html.Hr(),
                        html.Strong("Mass spectrum"),
                        html.Div(
                            dcc.Graph(
                                id="ms-spectrum",
                                figure=go.Figure(
                                    layout=go.Layout(
                                        title="Mass Spectrum",
                                        xaxis=dict(title="m/z", range=[0, 150]),
                                        yaxis=dict(title="Relative Intensity", range=[0, 1]),
                                        template="plotly_white"
                                    )
                                )
                            )
                        ),
                        html.Hr(style={"border-top": "5px solid #00305d"}),
                    ]
                )
            ]
        ),
        

        ############### Impressum ###############
        dbc.Container(
            [
                
                html.Div(
                    [
                        # button to expand/collapse
                        html.Div(
                            [
                                html.Button(
                                    "Legal notice", 
                                    id="impressum-button", 
                                    n_clicks=0, 
                                    style={'font-size': '12px'}
                                ),
                            ],
                            style={
                                "display": "flex",
                                "justify-content": "center",
                                "align-items": "center"
                            }
                        ),

                        html.Br(),
                        html.Br(),

                        # foldable area
                        dbc.Collapse(
                            html.Div(
                                dcc.Markdown(
                                    """
                                    The [Legal Notice of TU Dresden](https://tu-dresden.de/impressum) applies with the following amendments:

                                    RESPONSIBILITIES 

                                    If you have any questions regarding content, please contact:  
                                    Lina Müller  
                                    Technische Universität Dresden  
                                    DE – 01062 Dresden  
                                    Email: lina.mueller@tu-dresden.de  
                                    Tel.: +49 351 463-32616  
                                    
                                    Technical implementation:  
                                    Technische Universität Dresden  
                                    Professur für Lebensmittelkunde und Bedarfsgegenstände  
                                    Bergstraße 66, DE – 01062 Dresden  
                                    Lina Müller  
                                    Email: lina.mueller@tu-dresden.de   

                                    
                                    DATA PROTECTION DECLARATION

                                    TU Dresden processes personal data for the use of the public website. This personal data pertains to cookies only, which are used exclusively for providing this service. In particular, this means that this website uses no tracking cookies to record or analyze user movement and behavior on our website. 
                                    
                                    Legal basis  
                                    The legal basis for this is Art. 6 para. 1 letter f GDPR.
                                    
                                    Rights of data subjects   
                                    —	You have the right to obtain information from TU Dresden on the data processed concerning you and/or to request the correction of inaccurate data.  
                                    —	You have the right to erasure and restriction of processing as well as the right to object to the processing.  
                                    —	You can contact TU Dresden's Data Protection Officer at any time:
                                    
                                    Technische Universität Dresden  
                                    Data Protection Officer  
                                    DE - 01062 Dresden  
                                    Tel.: +49 351 463 32839  
                                    Fax : +49 351 463 39718  
                                    Email: informationssicherheit@tu-dresden.de  
                                    https://tu-dresden.de/informationssicherheit  

                                    —	You also have the right to appeal to the supervisory authority if you believe that the processing of data concerning your person does not comply with the law. The supervisory authority for data protection is:
                                    Saxon Data Protection and Transparency Officer:

                                    Dr. Juliane Hundert  
                                    Maternistraße 17  
                                    DE - 01067 Dresden  
                                    Email: post@sdtb.sachsen.de   
                                    Phone: + 49 (0) 35185471 101  
                                    www.datenschutz.sachsen.de   
                                    """,
                                    link_target="_blank"  # ensures that links open in a new tab
                                ),
                                style={
                                    'font-size': '12px',
                                    "padding": "10px", 
                                    "border": "1px solid #ddd", 
                                    "borderRadius": "3px"
                                },
                            ),
                            id="impressum-collapse",
                            is_open=False
                        ),
                        html.Br(),
                        html.Br(),
                        html.Br()
                    ]
                ),
            ]
        ),
    ]
)



###############################################################################################
# CALLBACKS
###############################################################################################

##### Download Bibliography
@app.callback(
    Output("download-bib", "data"),
    Input("btn-download-bib", "n_clicks"),
    prevent_initial_call=True
)
def download_bib(n_clicks):
    if not n_clicks:
        return dash.no_update

    return dcc.send_file("Bibliography.txt")

##### Master table 
@app.callback(
    Output("table-count", "children"),
    Input("table-master", "derived_virtual_data"),
)
def update_count(rows):
    if rows is None:
        return "0 enteries"

    return f"Currently {len(rows)} entries in the table:"

@app.callback(
    Output("selected-substances", "data"),
    Input("btn-add", "n_clicks"),
    Input("btn-remove", "n_clicks"),
    State("table-master", "derived_virtual_data"),
    State("table-master", "derived_virtual_selected_rows"),
    State("table-selected", "selected_rows"),
    State("selected-substances", "data"),
    prevent_initial_call=True
)
def update_store(n_add, n_remove, master_data, master_selected, selected_rows_table, store):
    if store is None:
        store = []

    ctx = callback_context
    if not ctx.triggered:
        return store
    
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "btn-add":
        
        if master_selected and master_data:
            cas_list = [master_data[i]["CAS-No."] for i in master_selected]
            sim_dict = get_sim_for_cas(cas_list)
            
            for i in master_selected:
                row = master_data[i].copy()
                sim_string = sim_dict.get(row["CAS-No."], "")

                if sim_string:
                    sim_sorted = ", ".join(
                        map(str,
                            sorted(
                                int(float(f.strip()))
                                for f in sim_string.split(",")
                            )
                        )
                    )
                else:
                    sim_sorted = ""
                row["Fragments"] = sim_sorted
                if not any(d["CAS-No."] == row["CAS-No."] for d in store):
                    store.append(row)

    elif button_id == "btn-remove":
        if selected_rows_table:
            store = [row for idx, row in enumerate(store) if idx not in selected_rows_table]

    return store

##### Output table
@app.callback(
    Output("table-selected", "data"),
    Input("selected-substances", "data"),
)
def update_selected_table(store):
    if not store:
        return []

    return store

@app.callback(
    Output("output-sim", "children"),
    Output("output-sim-1", "children"),
    Input("selected-substances", "data"),
)
def update_sim_lists(store):

    if not store:
        return None, None

    cas_list = [s["CAS-No."] for s in store]

    # all fragments
    sim_dict = get_sim_for_cas(cas_list)
    all_frags = sorted({
        int(float(f.strip()))
        for v in sim_dict.values()
        for f in v.split(",")
    })

    # fragments with intensity 1.0
    sim1_dict = get_sim_for_cas(cas_list, highest_intensity_only=True)

    top_frags = sorted({
        int(float(f.strip()))
        for v in sim1_dict.values()
        for f in v.split(",")
    })

    return (
        html.Div([
            f"The fragment masses of all selected substances are combined here: ",
            ", ".join(map(str, all_frags))
        ]),
        html.Div([
            f"The masses of the respective fragments with the highest intensity are combined here: ",
            ", ".join(map(str, top_frags))
        ])
    )

##### Download output table
@app.callback(
    Output("download-excel", "data"),
    Input("btn-download-csv", "n_clicks"),
    State("selected-substances", "data"),
    prevent_initial_call=True
)
def download_csv(n_clicks, store):

    if not store:
        return dash.no_update

    df = pd.DataFrame(store)

    return dcc.send_data_frame(
        df.to_csv,
        "output-table.csv",
        index=False
    )

@app.callback(
    Output("download-csv", "data"),
    Input("btn-download-excel", "n_clicks"),
    State("selected-substances", "data"),
    prevent_initial_call=True
)
def download_excel(n_clicks, store):

    if not store:
        return dash.no_update

    df = pd.DataFrame(store)

    return dcc.send_data_frame(
        df.to_excel,
        "output-table.xlsx",
        index=False
    )

##### Chromatogram
@app.callback(
    Output("chromatogram", "figure"),
    Output("phase-warning", "children"),
    Input("btn-analyse", "n_clicks"),
    State("selected-substances", "data"),
    State("dropdown-phase", "value"),
    prevent_initial_call=True
)
def run_analysis(n_clicks, store, phase):
    if not store:
        return show_chrom([]), None
    
    if not phase:
        return show_chrom([]), "Please select a phase!"
    
    ri_list = get_ri_list(store, phase)
    return show_chrom(ri_list), None

##### Mass spectrum
@app.callback(
    Output("ms-spectrum", "figure"),
    Input("chromatogram", "clickData"),
    State("dropdown-phase", "value"),
    State("selected-substances", "data"),
    prevent_initial_call=True
)
def display_mass_spectrum(clickData, phase, store):
    if not clickData or not store or not phase:
        return go.Figure()

    clicked_ri = clickData["points"][0]["x"]

    # Find the substance whose RI is closest to
    col = phase
    df_store = pd.DataFrame(store)
    df_store["RI_diff"] = abs(df_store[col] - clicked_ri)
    closest_row = df_store.loc[df_store["RI_diff"].idxmin()]

    cas = closest_row["CAS-No."]
    name = closest_row["Odorant"]

    fig_ms = show_ms(cas, title=f"Mass Spectrum {name}")
    return fig_ms

##### impressum

@app.callback(
    Output("impressum-collapse", "is_open"),
    Input("impressum-button", "n_clicks"),
    State("impressum-collapse", "is_open")
)
def impressum(n, is_open):
    if n:
        return not is_open
    return is_open

###############################################################################################

if __name__ == "__main__":
    app.run(debug=False)