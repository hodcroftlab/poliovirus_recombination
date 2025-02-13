import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import pandas as pd
import plotly.express as px
import argparse
import fastcluster
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from scipy.signal import find_peaks

# Define the command-line arguments
argparser = argparse.ArgumentParser(description="Plot pairwise p-distances between genomic regions.")
argparser.add_argument("--pdistances", nargs="+", help="Input file(s) (.csv) containing pairwise p-distances.")
argparser.add_argument("--regions", nargs="+", help="Genomic regions to compare.")
argparser.add_argument("--metadata", help="Metadata file containing additional information on each sequence.")
argparser.add_argument("--cutoff", type=float, help="Cutoff value for p-distance difference.")
args = argparser.parse_args()

# Load data files
pdistances = args.pdistances
pdistance_dict = {} # Prepare the dictionary with dataframes
for file in pdistances:
    serotype = file.split("/")[-1].replace("_pdistances.csv", "")
    pdistance_dict[serotype] = pd.read_csv(file)

# List of genomic regions (assumed to be the same across all datasets)
regions = args.regions

# Initialize the Dash app
app = dash.Dash(__name__)

# Layout of the app
app.layout = html.Div([
    html.H1("Interactive P-Distance Plots", style={"font-family": "Arial"}),

    # All three dropdowns in the same row
    html.Div([
        # Dropdown for serotype selection
        html.Div([
            html.Label("Select Serotype:", style={"font-family": "Arial", "font-size": "18px"}),
            dcc.Dropdown(
                id="serotype-dropdown",
                options=[{"label": serotype, "value": serotype} for serotype in pdistance_dict.keys()],
                value=list(pdistance_dict.keys())[0],
                style={"font-family": "Arial", "font-size": "16px"}
            )
        ], style={"display": "inline-block", "width": "32%", "padding-right": "10px", "vertical-align": "top"}),

        # Dropdown for region1 selection
        html.Div([
            html.Label("Select Genomic Region 1:", style={"font-family": "Arial", "font-size": "18px"}),
            dcc.Dropdown(
                id="region1-dropdown",
                options=[{"label": region, "value": region} for region in regions],
                value=regions[0],
                style={"font-family": "Arial", "font-size": "16px"}
            )
        ], style={"display": "inline-block", "width": "32%", "padding-right": "10px", "vertical-align": "top"}),

        # Dropdown for region2 selection
        html.Div([
            html.Label("Select Genomic Region 2:", style={"font-family": "Arial", "font-size": "18px"}),
            dcc.Dropdown(
                id="region2-dropdown",
                options=[{"label": region, "value": region} for region in regions],
                value=regions[1],
                style={"font-family": "Arial", "font-size": "16px"}
            )
        ], style={"display": "inline-block", "width": "32%", "vertical-align": "top"}),
        
        # Dropdown to choose coloring method (diff_abs or cluster)
        html.Div([
            html.Label("Color By:", style={"font-family": "Arial", "font-size": "18px"}),
            dcc.Dropdown(
                id="color-by-dropdown",
                options=[
                    {"label": "Abs. P-Distance Difference", "value": "diff_abs"},
                    {"label": "Cluster Assignment", "value": "cluster"}
                ],
                value="diff_abs",  # Default color by diff_abs
                style={"font-family": "Arial", "font-size": "16px"}
            )
        ], style={"display": "inline-block", "width": "50%", "vertical-align": "top"})
    ]),

    # First plot (scatter plot)
    dcc.Graph(id="pdist-plot"),
    
    # Second plot (histogram)
    dcc.Graph(id="diff-distribution-plot")
], style={"display": "block", "margin-left": "auto", "margin-right": "auto", "width": "700px"})


# Callback to update region2 options based on region1 selection
@app.callback(
    Output("region2-dropdown", "options"),
    [Input("region1-dropdown", "value")]
)
def update_region2_options(selected_region1):
    # Filter out the selected region1 from the available options for region2
    return [{"label": region, "value": region} for region in regions if region != selected_region1]

# Callback to update the plot based on dropdown selections
@app.callback(
    [Output("pdist-plot", "figure"), Output("diff-distribution-plot", "figure")],
    [Input("serotype-dropdown", "value"),
     Input("region1-dropdown", "value"),
     Input("region2-dropdown", "value"),
     Input("color-by-dropdown", "value")]
)
def update_plot(selected_serotype, region1, region2, color_by):
    # Select the dataframe for the chosen serotype
    df = pdistance_dict[selected_serotype]
    
    # Create a new dataframe for the selected regions
    df_subset = df[["seq1", "seq2", region1, region2]].copy()
    
    # Remove rows where with missing values
    df_subset = df_subset.dropna()
    
    # Calculate difference between region1 and region2
    if region1 == "5UTR":
        df_subset["diff"] = (df_subset[region1] - df_subset[region2] * 0.8).round(4)
    elif region2 == "5UTR":
        df_subset["diff"] = (df_subset[region1] * 0.8 - df_subset[region2]).round(4)
    else:
        df_subset["diff"] = (df_subset[region1] - df_subset[region2]).round(4)
                
    df_subset["diff_abs"] = df_subset["diff"].abs()
    
    # Apply agglomerative clustering to the p-distance plot
    X = np.array(df_subset[[region1, region2]], dtype=np.float64)
    clustering = AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=0.01)
    df_subset["classification"] = clustering.fit_predict(X)
    
    color_variable = "diff_abs" if color_by == "diff_abs" else "classification"
    
    
    # Determine max value for the axes
    max_val = max(df_subset[region1].max(), df[region2].max())
    
    # Create the plot
    fig1 = px.scatter(
        df_subset,
        x=region1,
        y=region2,
        color=color_variable,
        hover_data=["seq1", "seq2", region1, region2, "diff"],
        title=f"{selected_serotype}: {region1} vs. {region2}<br><sub>Number of unique sequences: {df_subset['seq1'].nunique()}</sub>",
        labels={region1: f"{region1} p-distance", region2: f"{region2} p-distance"},
        color_continuous_scale="Viridis" if color_by == "diff_abs" else "Rainbow"
    )
    
    # Add lines for reference (diagonal) and cut-offs
    cutoff = args.cutoff
    
    # Adjust the line slope conditionally if 5'UTR is chosen as region1 or region2
    if region1 == "5UTR":
        slope = 1/0.8  # Higher slope when 5'UTR is on the y-axis

        fig1.add_shape(
            type="line",
            x0=0, y0=0,
            x1=1, y1=slope, 
            line=dict(color="black", dash="dash", width=1)
        )
        
        fig1.add_shape(
            type="line",
            x0=-cutoff, y0=0,
            x1=1-cutoff, y1=slope,
            line=dict(color="red", dash="dot", width=0.5)
        )
        
        fig1.add_shape(
            type="line",
            x0=cutoff, y0=0,
            x1=1+cutoff, y1=slope,
            line=dict(color="red", dash="dot", width=0.5)
        )
        
    elif region2 == "5UTR":
        slope = 0.8  # Lower slope when 5'UTR is on the x-axis
        
            # Add the diagonal reference line with appropriate slope
        fig1.add_shape(
            type="line",
            x0=0, y0=0,
            x1=1, y1=slope, 
            line=dict(color="black", dash="dash", width=1)
        )
        
        fig1.add_shape(
            type="line",
            x0=0, y0=cutoff,
            x1=1, y1=slope+cutoff,
            line=dict(color="red", dash="dot", width=0.5)
        )
        
        fig1.add_shape(
            type="line",
            x0=0, y0=-cutoff,
            x1=1, y1=slope-cutoff,
            line=dict(color="red", dash="dot", width=0.5)
        )
        
    else: 
        slope = 1 # Default slope
        
            # Add the diagonal reference line with appropriate slope
        fig1.add_shape(
            type="line",
            x0=0, y0=0,
            x1=1, y1=slope, 
            line=dict(color="black", dash="dash", width=1)
        )
        
        fig1.add_shape(
            type="line",
            x0=0, y0=cutoff,
            x1=1, y1=slope+cutoff,
            line=dict(color="red", dash="dot", width=0.5)
        )
        
        fig1.add_shape(
            type="line",
            x0=cutoff, y0=0,
            x1=1, y1=slope-cutoff,
            line=dict(color="red", dash="dot", width=0.5)
        )
    
    # Update the layout to make the plot quadratic, black axes, no grid
    fig1.update_layout(
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            color="black",
            range=[-0.01, 1.01 * max_val],
            showline=True,
            linecolor="black",
            ticks="outside",
            tickwidth=2
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            color="black",
            range=[-0.01, 1.01 * max_val],
            showline=True,
            linecolor="black",
            ticks="outside",
            tickwidth=2
        ),
        plot_bgcolor="white",
        autosize=False,
        width=700,
        height=700
    )
    
    # Calculate mean and standard deviation of the "diff" values
    mean_diff = df_subset["diff"].mean()
    std_diff = df_subset["diff"].std()
    
    
    # Create the second plot: a histogram of diff values
    fig2 = px.histogram(
        df_subset,
        x="diff",
        nbins=100, 
        marginal="violin",  # Add small violin plot to show distribution
        title=f"{selected_serotype}: P-Distance Differences ({region1} vs. {region2})",
        labels={"diff": f"Difference in P-Distances ({region1} - {region2})"}
    )
    
    # Add vertical lines for midpoint, mean, standard deviation, and cutoff
    fig2.add_shape(
        type="line",
        x0=0, y0=0, x1=0, y1=1,
        line=dict(color="black", width=1),
        xref="x", yref="paper",
        name="Midpoint"
    )
    
    fig2.add_shape(
        type="line",
        x0=mean_diff, y0=0, x1=mean_diff, y1=1,
        line=dict(color="grey", width=1),
        xref="x", yref="paper", 
        name="Mean"
    )
    fig2.add_shape(
        type="line",
        x0=mean_diff - std_diff, y0=0, x1=mean_diff - std_diff, y1=1,
        line=dict(color="grey", width=1, dash="dot"),
        xref="x", yref="paper",
        name="-1 SD"
    )
    fig2.add_shape(
        type="line",
        x0=mean_diff + std_diff, y0=0, x1=mean_diff + std_diff, y1=1,
        line=dict(color="grey", width=1, dash="dot"),
        xref="x", yref="paper",
        name="+1 SD"
    )
    
    fig2.add_shape(
        type="line",
        x0=cutoff, y0=0, x1=cutoff, y1=0.86,
        line=dict(color="red", width=1, dash="dot"),
        xref="x", yref="paper",
        name="Cutoff"
    )
    
    fig2.add_shape(
        type="line",
        x0=-cutoff, y0=0, x1=-cutoff, y1=0.86,
        line=dict(color="red", width=1, dash="dot"),
        xref="x", yref="paper",
        name="Cutoff"
    )
    
    # Add text labels above the lines for mean, ±1 SD, and cutoff (rotated by 90 degrees)
    fig2.add_annotation(
        x=mean_diff, y=1.05,  # Position above the plot
        text="Mean",
        showarrow=False,
        font=dict(color="grey", size=10),
        xref="x", yref="paper"
    )
    fig2.add_annotation(
        x=mean_diff - std_diff, y=1.05,
        text="-1 SD",
        showarrow=False,
        font=dict(color="grey", size=10),
        xref="x", yref="paper"
    )
    fig2.add_annotation(
        x=mean_diff + std_diff, y=1.05,
        text="+1 SD",
        showarrow=False,
        font=dict(color="grey", size=10),
        xref="x", yref="paper"
    )
    fig2.add_annotation(
        x=cutoff-0.005, y=0.8,
        text="Cutoff",
        textangle=-90,  # Rotate the text by 90 degrees
        showarrow=False,
        font=dict(color="red", size=10),
        xref="x", yref="paper"
    )
    fig2.add_annotation(
        x=-cutoff-0.005, y=0.8,
        text="Cutoff",
        textangle=-90,  # Rotate the text by 90 degrees
        showarrow=False,
        font=dict(color="red", size=10),
        xref="x", yref="paper"
    )
    
    fig2.update_layout(
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            color="black",
            showline=True,
            linecolor="black",
            ticks="outside",
            tickwidth=2
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            color="black",
            showline=True,
            linecolor="black",
            ticks="outside",
            tickwidth=2
        ),
        plot_bgcolor="white",
        width=700,
        height=600
    )

    return fig1, fig2

# Export to HTML
if __name__ == "__main__":
    # Run the app locally
    app.run_server(debug=False, port=8050)
