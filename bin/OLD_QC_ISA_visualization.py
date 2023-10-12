#!/bin/env python

# %%
# std imports
import argparse
from typing import Optional
import math
from io import BytesIO
from pathlib import Path
import zipfile

# 3rd party imports
import pandas as pd
import numpy as np
from sqlalchemy import create_engine, text, bindparam
import plotly
import plotly.express as px
import plotly.io as pio
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

pio.renderers.default = "png"

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-csv_file", help="CSV file with QC data.")
    parser.add_argument("-output", help="Output folder for the plots as json files.")
    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_setup()

    #############
    # parameters

    ### use some ISA runs for testing
    csv_file = args.csv_file.split(",")
    #csv_file = "temp/example_nf_out.csv"

    ### folder to save the plots as json files
    output_path = args.output    # "graphics/"
    #output_path = "graphics/"

    ##########
    ISA_df = pd.read_csv(csv_file)

  

    # Figure 1: Number of proteins, protein groups and unfiltered protein groups
    #### TODO: has to be tested, when PIA output is finished!
    df_pl1 = ISA_df[["filename", "nrProteins", "number-filtered-protein-groups", "nrProteingroups_unfiltered"]]
    df_pl1_long = df_pl1.melt(id_vars = ["filename"])
    # TODO: evtl. Spalten umbenennen damit es schoener aussieht?

    fig1 = px.bar(df_pl1_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of proteins, protein groups and unfiltered protein groups")
    fig1.update_yaxes(exponentformat="none") 

    with open(output_path +"/isafig1_barplot_proteins_proteingroups.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig1))


    # Figure 2: Number of peptides
    #### TODO: has to be tested, when PIA output is finished!
    df_pl2 = ISA_df[["filename", "number-filtered-peptides"]]

    fig2 = px.bar(df_pl2, x="filename", y="number-filtered-peptides", #color="variable", barmode = "group", 
                title = "Number of peptides")
    fig2.update_yaxes(exponentformat="none") 

    with open(output_path +"/isafig2_barplot_peptides.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig2))

    
    # Figure 3: Number of PSMs
    #### TODO: has to be tested, when PIA output is finished!
    df_pl3 = ISA_df[["filename", "number-filtered-psms"]]

    fig3 = px.bar(df_pl3, x="filename", y="number-filtered-psms", #color="variable", barmode = "group", 
                title = "Number of PSMs")
    fig3.update_yaxes(exponentformat="none") 

    with open(output_path +"/isafig3_barplot_PSMs.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig3))


    # Figure 4: 
    ## TODO: add ion map plot
    ### Time vs. Lockmass deviation
    ### pressure (over time)
    ### ion injection time (over retention time)


