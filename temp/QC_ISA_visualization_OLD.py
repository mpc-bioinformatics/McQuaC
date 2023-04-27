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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-raw_files", help="List of raw files to plot from database (comma-separated).")
    parser.add_argument("-output", help="Output folder for the plots as json files.")
    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_setup()

    #############
    # parameters

    ### use some ISA runs for testing
    raw_files = args.raw_files.split(",")
    #raw_files = ["FLI18416std", "FLI18414std", "EXII01692std", "QEXI39066std"]

    ### folder to save the plots as json files
    output_path = args.output    # "graphics/"
    #output_path = "graphics/"


    ##########

    engine = create_engine("mysql+mariadbconnector://mpcqc:quality@mpc-qc/mpcqc")

    file_df: Optional[pd.DataFrame] = None
    with engine.connect() as conn:
        statement = text("SELECT * FROM files WHERE filename IN :raw_files")
        statement = statement.bindparams(
            bindparam("raw_files", tuple(raw_files), expanding=True)
        )
        query = conn.execute(statement)
        file_df = pd.DataFrame(query.fetchall(), columns=query.keys())

    ids = file_df["id"]

    ISA_df: Optional[pd.DataFrame] = None
    with engine.connect() as conn:
        statement = text("SELECT * FROM isa_data WHERE fileId IN :ids")
        statement = statement.bindparams(
            bindparam("ids", tuple(ids), expanding=True)
        )
        print(statement.compile(compile_kwargs={"literal_binds": True}))
        query = conn.execute(statement)
        ISA_df = pd.DataFrame(query.fetchall(), columns=query.keys())
        
    #print(ISA_df)

    # %%
    # Figure 1: Number of proteins, protein groups and unfiltered protein groups

    df_pl1 = ISA_df[["fileId", "run_name", "nrProteins", "nrProteingroups", "nrProteingroups_unfiltered"]]
    df_pl1_long = df_pl1.melt(id_vars = ["fileId", "run_name"])
    #print(df_pl1_long)
    #df_pl1_long = df_pl1_long.reset_index(level=["level"])

    fig1 = px.bar(df_pl1_long, x="run_name", y="value", color="variable", barmode = "group", 
                title = "Number of proteins, protein groups and unfiltered protein groups")
    fig1.update_yaxes(exponentformat="none") 
    #fig1.show()

    with open(output_path +"/isafig1_barplot_proteins_proteingroups.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig1))

    # %%
    # Figure 2: Number of peptides

    df_pl2 = ISA_df[["fileId", "run_name", "nrPeptides"]]
    #df_pl1_long = df_pl1.melt(id_vars = ["fileId", "run_name"])
    #print(df_pl1_long)
    #df_pl1_long = df_pl1_long.reset_index(level=["level"])

    fig2 = px.bar(df_pl2, x="run_name", y="nrPeptides", #color="variable", barmode = "group", 
                title = "Number of peptides")
    fig2.update_yaxes(exponentformat="none") 
    #fig2.show()

    with open(output_path +"/isafig2_barplot_peptides.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig2))

    # %%
    # Figure 3: Number of PSMs

    df_pl3 = ISA_df[["fileId", "run_name", "nrPSMs"]]
    #df_pl1_long = df_pl1.melt(id_vars = ["fileId", "run_name"])
    #print(df_pl1_long)
    #df_pl1_long = df_pl1_long.reset_index(level=["level"])

    fig3 = px.bar(df_pl3, x="run_name", y="nrPSMs", #color="variable", barmode = "group", 
                title = "Number of PSMs")
    fig3.update_yaxes(exponentformat="none") 
    #fig3.show()

    with open(output_path +"/isafig3_barplot_PSMs.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig3))


