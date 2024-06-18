#!/usr/bin/env python

# %%
import re
# std imports
import argparse
#from typing import Optional
#import math
#from io import BytesIO
#from pathlib import Path
#import zipfile
import ast
import os
import base64
import pickle
import zlib
import datetime

# 3rd party imports
import pandas as pd
import numpy as np
#from sqlalchemy import create_engine, text, bindparam
import plotly
import plotly.express as px
import plotly.graph_objects as go
#import plotly.offline as pyo
import plotly.io as pio
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

pio.renderers.default = "png"


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-csv_file", help="CSV file with QC data.", default=None)
    parser.add_argument("-group", help="List of the experimental group (comma-separated).", default=None)
    parser.add_argument("-output", help="Output folder for the plots as json files.", default = "graphics")
    parser.add_argument("-tic_overlay_offset", help = "Offset for TIC overlay plots", default = 0)
    parser.add_argument("-fig_show", help = "Show figures, e.g. for debugging?", default = False, action = "store_true")
    parser.add_argument("-isa", help = "Is this ISA QC?", default = False, action = "store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()


####################################################################################################
    # parameters

    ### use HCC dataset for testing
    csv_file = args.csv_file#.split(",")

    ### grouping
    ### If use_groups = False, PCA plots are coloured by timestamp. If True, PCA plots are coloured by group.
    if (args.group is None):
        use_group = False
    else:
        use_group = True
        group = np.array(args.group.split(","))
    
    ### folder to save the plots as json files
    output_path = args.output    

    fig_show = args.fig_show
    isa = args.isa

    #csv_file = "temp/quality_control_20230920.csv"

####################################################################################################
    # read csv file
    df = pd.read_csv(csv_file)
    nr_rawfiles = len(df)

    # add filenames if not present
    def add_filename_column(dataframe, column_name):
        if column_name not in dataframe.columns:
            filenames = [f"file{i}" for i in range(len(dataframe))]
            dataframe.insert(0, column_name, filenames)
        return dataframe
    df = add_filename_column(df, "filename")
    
    # sort data by name run name to have the right order for the plots
    df = df.sort_values(by = "filename")
    df.reset_index(drop = True, inplace = True)
    

####################################################################################################
    # short function to extract infos from compressed columns
    def unbase64_uncomp_unpickle(x: bytes) -> list:
        unb64 = base64.b64decode(x)
        uncomp = zlib.decompress(unb64)
        return pickle.loads(uncomp)


    ############################################################################################
    # Table 0: table with overview over QC measures. 
    feature_list = ["file_and_analysis_timestamp",
                    "total_num_ms1",
                    "total_num_ms2",
                    "number_filtered_psms",
                    "number_filtered_peptides", 
                    "number_proteins", 
                    "total_num_ident_features",
                    "total_num_features",
                    "RT_duration", 
                    "RT_TIC_Q_000-025",
                    "RT_TIC_Q_025-050",
                    "RT_TIC_Q_050-075",
                    "RT_TIC_Q_075-100",
                    "RT_MS1_Q_000-025",
                    "RT_MS1_Q_025-050",
                    "RT_MS1_Q_050-075",
                    "RT_MS1_Q_075-100",
                    "RT_MS2_Q_000-025",
                    "RT_MS2_Q_025-050", 
                    "RT_MS2_Q_050-075",
                    "RT_MS2_Q_075-100", 
                    "MS1-TIC-Change-Q2",
                    "MS1-TIC-Change-Q3", 
                    "MS1-TIC-Change-Q4",
                    "MS1-TIC-Q2",
                    "MS1-TIC-Q3",
                    "MS1-TIC-Q4", 
                    "MS1_Freq_Max",
                    "MS1_Density_Q1",
                    "MS1_Density_Q2", 
                    "MS1_Density_Q3",
                    "MS2_Freq_Max",
                    "MS2_Density_Q1",
                    "MS2_Density_Q2", 
                    "MS2_Density_Q3",
                    "MS2_PrecZ_1", 
                    "MS2_PrecZ_2", 
                    "MS2_PrecZ_3",
                    "MS2_PrecZ_4",
                    "MS2_PrecZ_5",
                    "MS2_PrecZ_more", 
                    "accumulated-MS1_TIC", 
                    "accumulated-MS2_TIC",
                    "num_features_charge_1",
                    "num_features_charge_2",
                    "num_features_charge_3",
                    "num_features_charge_4",
                    "num_features_charge_5", 
                    "psm_charge1",
                    "psm_charge2", 
                    "psm_charge3", 
                    "psm_charge4", 
                    "psm_charge5",
                    "psm_charge_more", 
                    "psm_missed_0",
                    "psm_missed_1",
                    "psm_missed_2",
                    "psm_missed_3",
                    "psm_missed_more"   
                    # TODO we should add then dynamically. Iow, we could thest if it contains binary data and if not include it into the list
                    # TODO These below are Thermo TUNE Information which is one dimensional
                    # "THEREMO_TUNE_ '''Ion Transfer Tube Temperature (+ or +-)'''"
                    # "THEREMO_TUNE_ '''Ion Transfer Tube Temperature (-)'''"
                    # "THEREMO_TUNE_ '''Vaporizer Temp. (+ or +-)'''"
                    # "THEREMO_TUNE_ '''Vaporizer Temp. (-)'''"
                ]
    
    x = [datetime.datetime.utcfromtimestamp(x) for x in df["timestamp"]] # convert timestamp to datetime
    
    ### for now, just check if the items in featurelist are really in the data frame and if not, skip them
    valid_features = [feature for feature in feature_list if feature in df.columns]
    df_table0 = df[valid_features]
    df_table0 = df_table0.loc[:,:].copy()
    df_table0.loc[:,"timestamp"] = x
    
    ### if file ISAs are analyzed, add the spike-ins to the table
    
    if isa:
        spikes = []
        for index in df.index:
            spike_tmp = unbase64_uncomp_unpickle(df["MPCSPIKEINS_____pickle_zlib"].iloc[index])
            spikes.append(spike_tmp)
        spikes_df = pd.DataFrame(spikes)
        df_table0 = pd.concat([df_table0, spikes_df], axis = 1)
        
    df_table0.to_csv(output_path + "/table0_summary.csv")

################################################################################################
    # Figure 1: Barplot for total number of MS1 and MS2 spectra
    df_pl1 = df[["filename", "total_num_ms1", "total_num_ms2"]]
    df_pl1_long = df_pl1.melt(id_vars = ["filename"])
    fig1 = px.bar(df_pl1_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Total number of MS1 and MS2 spectra")
    fig1.update_yaxes(exponentformat="none") 
    fig1.update_xaxes(tickangle=-90)
    if fig_show: 
        fig1.show()
    with open(output_path +"/fig1_barplot_MS1_MS2.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig1))
    #pyo.plot(fig1, filename = output_path +"/fig1_barplot_MS1_MS2.html")
    fig1.write_html(file = output_path +"/fig1_barplot_MS1_MS2.html", auto_open = False)


################################################################################################
    # Figure 2: Barplot for number of PSMs, peptides, proteins
    
    if ("number_filtered_psms" in df.columns):  # e.g. for DIA, no identification is done so those columns are missing
        df_pl2 = df[["filename", "number_filtered_psms", "number_filtered_peptides", "number_proteins", "number_ungrouped_proteins"]]
        df_pl2_long = df_pl2.melt(id_vars = ["filename"])
        fig2 = px.bar(df_pl2_long, x="filename", y="value", color="variable", barmode = "group", 
                    title = "Number of filtered PSMs, filtered peptides, filtered protein groups and accessions")
        fig2.update_yaxes(exponentformat="none") 
        fig2.update_xaxes(tickangle=-90)
    else: 
        fig2 = go.Figure()
        fig2.add_annotation(
            x=0.5,
            y=0.5,
            text="Columns are missing, no plot created!",
            showarrow=False,
            font=dict(size=14)
        )
        fig2.update_layout(
            width=600,
            height=400,
            title="Empty Plot"
        )
    if fig_show: 
        fig2.show()
    with open(output_path +"/fig2_barplot_PSMs_peptides_proteins.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig2))
    #pyo.plot(fig2, filename = output_path +"/fig2_barplot_PSMs_peptides_proteins.html")
    fig2.write_html(file = output_path +"/fig2_barplot_PSMs_peptides_proteins.html", auto_open = False)
    
################################################################################################
    # Figure 3: Barplot for features and identified features
    df_pl3 = df[["filename", "total_num_features", "total_num_ident_features"]]
    df_pl3_long = df_pl3.melt(id_vars = ["filename"])
    fig3 = px.bar(df_pl3_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of features and identified features")
    fig3.update_yaxes(exponentformat="none") 
    fig3.update_xaxes(tickangle=-90)
    if fig_show: 
        fig3.show()
    with open(output_path +"/fig3_barplot_features.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig3))
    #pyo.plot(fig3, filename = output_path +"/fig3_barplot_features.html")
    fig3.write_html(file = output_path +"/fig3_barplot_features.html", auto_open = False)


####################################################################################################
    ## Figure 4: TIC Overlay as Lineplot
    offset_ratio = float(args.tic_overlay_offset)

    tic_df = []
    for index in df.index:
        tic = dict(TIC=ast.literal_eval(df["ms1_tic_array"].iloc[index]),
            RT=ast.literal_eval(df["ms1_rt_array"].iloc[index]),
            filename=[df["filename"].iloc[index]]*len(ast.literal_eval(df["ms1_tic_array"].iloc[index])))
        tic = pd.DataFrame(tic)
        tic_df.append(tic)
    tic_df2 = pd.concat(tic_df)

    offset = max(tic_df2["TIC"])


    offset_tmp = 0
    tic_df3 = []
    for tmp in tic_df:
        tmp["TIC"] = tmp["TIC"] + offset_tmp
        tic_df3.append(tmp)
        offset_tmp = offset_tmp + offset_ratio*offset
    tic_df = pd.concat(tic_df3)

    fig4 = px.line(tic_df, x="RT", y="TIC", color = "filename", title = "TIC overlay")
    fig4.update_traces(line=dict(width=0.5))
    fig4.update_yaxes(exponentformat="E") 
    if fig_show:
        fig4.show()
    with open(output_path +"/fig4_TIC_overlay.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig4))
    #pyo.plot(fig4, filename = output_path +"/fig4_TIC_overlay.html")
    fig4.write_html(file = output_path +"/fig4_TIC_overlay.html", auto_open = False)


#################################################################################################
    # Figure 5: Barplot TIC quartiles
    df_pl5 = df[["filename", 'RT_TIC_Q_000-025', 'RT_TIC_Q_025-050', 'RT_TIC_Q_050-075', 'RT_TIC_Q_075-100']]
    df_pl5_long = df_pl5.melt(id_vars = ["filename"])
    fig5 = px.bar(df_pl5_long, x = "filename", y = "value", color = "variable", title = "Quartiles of TIC over retention time")
    fig5.update_xaxes(tickangle=-90)
    if fig_show:
        fig5.show()
    with open(output_path +"/fig5_barplot_TIC_quartiles.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig5))
    #pyo.plot(fig5, filename = output_path +"/fig5_barplot_TIC_quartiles.html")
    fig5.write_html(file = output_path +"/fig5_barplot_TIC_quartiles.html", auto_open = False)

################################################################################################
    # Figure 6: Barplot MS1 TIC quartiles
    df_pl6 = df[["filename", 'RT_MS1_Q_000-025', 'RT_MS1_Q_025-050', 'RT_MS1_Q_050-075', 'RT_MS1_Q_075-100']]
    df_pl6_long = df_pl6.melt(id_vars = ["filename"])
    fig6 = px.bar(df_pl6_long, x="filename", y="value", color="variable", title = "Quartiles of MS1 over retention time")
    fig6.update_xaxes(tickangle=-90)
    if fig_show: 
        fig6.show()
    with open(output_path +"/fig6_barplot_MS1_TIC_quartiles.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig6))
    #pyo.plot(fig6, filename = output_path +"/fig6_barplot_MS1_TIC_quartiles.html")   
    fig6.write_html(file = output_path +"/fig6_barplot_MS1_TIC_quartiles.html", auto_open = False)

################################################################################################
    # Figure 7: Barplot MS2 TIC quartiles
    df_pl7 = df[["filename", 'RT_MS2_Q_000-025', 'RT_MS2_Q_025-050', 'RT_MS2_Q_050-075', 'RT_MS2_Q_075-100']]
    df_pl7_long = df_pl7.melt(id_vars = ["filename"])
    fig7 = px.bar(df_pl7_long, x="filename", y="value", color="variable", title = "Quartiles of MS2 over retention time")
    fig7.update_xaxes(tickangle=-90)
    if fig_show:
        fig7.show()
    with open(output_path +"/fig7_barplot_MS2_TIC_quartiles.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig7))
    #pyo.plot(fig7, filename = output_path +"/fig7_barplot_MS2_TIC_quartiles.html")
    fig7.write_html(file = output_path +"/fig7_barplot_MS2_TIC_quartiles.html", auto_open = False)

################################################################################################
    # Figure 8: Precursor charge states
    df_pl8 = df[["filename", 'MS2_PrecZ_1', 'MS2_PrecZ_2', 'MS2_PrecZ_3', 'MS2_PrecZ_4', 'MS2_PrecZ_5', 'MS2_PrecZ_more', 'MS2_PrecZ_Unknown']]
    df_pl8_long = df_pl8.melt(id_vars = ["filename"])
    fig8 = px.bar(df_pl8_long, x="filename", y="value", color="variable", title = "Charge states of precursors")
    fig8.update_xaxes(tickangle=-90)
    if fig_show:
        fig8.show()
    with open(output_path +"/fig8_barplot_precursor_chargestate.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig8))
    #pyo.plot(fig8, filename = output_path +"/fig8_barplot__precursor_chargestate.html")
    fig8.write_html(file = output_path +"/fig8_barplot__precursor_chargestate.html", auto_open = False)

################################################################################################
    # Figure 9: PSM charge states (of identified spectra)
    #### TODO: relative statt absolute Zahlen!
    
    if ("psm_charge1" in df.columns):
        df_pl9 = df[["filename", 'psm_charge1', 'psm_charge2', 'psm_charge3', 'psm_charge4', 'psm_charge5', 'psm_charge_more']]
        df_pl9_long = df_pl9.melt(id_vars = ["filename"])
        fig9 = px.bar(df_pl9_long, x="filename", y="value", color="variable", title = "Charge states of PSMs")
        fig9.update_xaxes(tickangle=-90)
    else: 
        fig9 = go.Figure()
        fig9.add_annotation(
            x=0.5,
            y=0.5,
            text="Columns are missing, no plot created!",
            showarrow=False,
            font=dict(size=14)
        )
        fig9.update_layout(
            width=600,
            height=400,
            title="Empty Plot"
        )
    if fig_show:
        fig9.show()
    with open(output_path +"/fig9_barplot_PSM_chargestate.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig9))
    #pyo.plot(fig9, filename = output_path +"/fig9_barplot_PSM_chargestate.html")
    fig9.write_html(file = output_path +"/fig9_barplot_PSM_chargestate.html", auto_open = False)

################################################################################################
    # Figure 10: Missed cleavages of PSMs
    ### TODO: relative statt absolute Zahlen!
    
    if ("psm_missed_0" in df.columns):
        df_pl10 = df[["filename", 'psm_missed_0', 'psm_missed_1', 'psm_missed_2', 'psm_missed_3', 'psm_missed_more']]
        df_pl10_long = df_pl10.melt(id_vars = ["filename"])
        fig10 = px.bar(df_pl10_long, x="filename", y="value", color="variable", title = "Number of missed cleavages for PSMs")
        fig10.update_xaxes(tickangle=-90)
    else: 
        fig10 = go.Figure()
        fig10.add_annotation(
            x=0.5,
            y=0.5,
            text="Columns are missing, no plot created!",
            showarrow=False,
            font=dict(size=14)
        )
        fig10.update_layout(
            width=600,
            height=400,
            title="Empty Plot"
        )  
    
    if fig_show:
        fig10.show()
    with open(output_path +"/fig10_barplot_PSM_missedcleavages.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig10))
    #pyo.plot(fig10, filename = output_path +"/fig10_barplot_PSM_missedcleavages.html")
    fig10.write_html(file = output_path +"/fig10_barplot_PSM_missedcleavages.html", auto_open = False)


#################################################################################################
    #### Preparations for PCA
    # If no groups are given, the points in the PCA are coloured by the timestamp.
    # The oldest raw file is coloured blue, the newest red.
    # Depending on the timestamp, the other ones are coloured on a gradient in between these two colours.

    timestamps = df[["timestamp"]].values.flatten().tolist()
    mintime = min(timestamps)
    maxtime = max(timestamps)

    ## Calculate scaling of timestamps to colour the points in the PCA plot (percentage between min and max time):
    t_scaled = []
    for t in timestamps:
        if (mintime == maxtime):
            t_scaled_tmp = 1
        else:
            t_scaled_tmp = (t - mintime)/(maxtime-mintime)*100 
        t_scaled.append(t_scaled_tmp)



    # Fig 11 PCA on all data
    if nr_rawfiles > 1:
        feature_list = ["RT_duration", 
        "total_num_ms1",
        "total_num_ms2", 
        "RT_TIC_Q_000-025",
        "RT_TIC_Q_025-050",
        "RT_TIC_Q_050-075",
        "RT_TIC_Q_075-100",
        "RT_MS1_Q_000-025",
        "RT_MS1_Q_025-050",
        "RT_MS1_Q_050-075",
        "RT_MS1_Q_075-100",
        "RT_MS2_Q_000-025",
        "RT_MS2_Q_025-050", 
        "RT_MS2_Q_050-075",
        "RT_MS2_Q_075-100", 
        "MS1-TIC-Change-Q2",
        "MS1-TIC-Change-Q3", 
        "MS1-TIC-Change-Q4",
        "MS1-TIC-Q2",
        "MS1-TIC-Q3",
        "MS1-TIC-Q4", 
        "MS1_Freq_Max",
        "MS1_Density_Q1",
        "MS1_Density_Q2", 
        "MS1_Density_Q3",
        "MS2_Freq_Max",
        "MS2_Density_Q1",
        "MS2_Density_Q2", 
        "MS2_Density_Q3",
        "MS2_PrecZ_1", 
        "MS2_PrecZ_2", 
        "MS2_PrecZ_3",
        "MS2_PrecZ_4",
        "MS2_PrecZ_5",
        "MS2_PrecZ_more", 
        "accumulated-MS1_TIC", 
        "accumulated-MS2_TIC",
        "total_num_ident_features",
        "num_features_charge_1",
        "num_features_charge_2",
        "num_features_charge_3",
        "num_features_charge_4",
        "num_features_charge_5", 
        "psm_charge1",
        "psm_charge2", 
        "psm_charge3", 
        "psm_charge4", 
        "psm_charge5",
        "psm_charge_more"
        ]

        ### for now, just check if the items in featurelist are really in the data frame and if not, skip them
        valid_features = [feature for feature in feature_list if feature in df.columns]
        
        df_pl11 = df[valid_features]
        df_pl11 = df_pl11.fillna(value = 0) # impute missig values by 0

        df_pl11_norm = pd.DataFrame(StandardScaler().fit_transform(df_pl11)) 

        #perform PCA
        pca = PCA(n_components=2)

        principalComponents = pca.fit_transform(df_pl11_norm)

        col = range(1,(principalComponents.shape[1]+1))
        col = ['pca'+ str(y) for y in col]
        principalDf = pd.DataFrame(data = principalComponents, columns = col)

        principalDf["t_scaled"] = t_scaled
        principalDf["raw_file"] = df["filename"]

        # Explained variance in axis labels
        pca_var = pca.explained_variance_ratio_
        pca_var1 = round(pca_var[0]*100, 1)
        pca_var2 = round(pca_var[1]*100, 1)

        label_x = "PC1 (" + str(pca_var1) + "%)"
        label_y = "PC2 (" + str(pca_var2) + "%)"

        if use_group:
            principalDf["group"] = group
            fig11 = px.scatter(principalDf, x = "pca1", y = "pca2", color = "group",
                        color_continuous_scale="bluered", title = "PCA on raw data",
                        hover_name="raw_file", hover_data=["pca1", "pca2"],
                            labels={
                            "pca1": label_x,
                            "pca2": label_y,
                            "group": "Group"
                        })
        else:
            principalDf["t_scaled"] = t_scaled
            fig11 = px.scatter(principalDf, x = "pca1", y = "pca2", color = "t_scaled",
                        color_continuous_scale="bluered", title = "PCA on all data",
                        hover_name="raw_file", hover_data=["pca1", "pca2"],
                            labels={
                            "pca1": label_x,
                            "pca2": label_y,
                            "t_scaled": "timestamp"
                        })
            
        # Table and plot with feature loadings (weights of the variables in the PCA)
        loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_features)
        loadings["length"] = np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2)
        loadings["variable"] = loadings.index
        fig11_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (all data)", 
            hover_name="variable", hover_data=["PC1", "PC2"],)
    else:
        fig11 = go.Figure()
        fig11.add_annotation(
            x=0.5,
            y=0.5,
            text="PCA cannot be computed using only one sample!",
            showarrow=False,
            font=dict(size=14)
        )
        fig11.update_layout(
            width=600,
            height=400,
            title="Empty Plot"
        )
        fig11_loadings = fig11
    if fig_show:
        fig11.show()
    with open(output_path +"/fig11a_PCA_all.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11))
    #pyo.plot(fig11, filename = output_path +"/fig11a_PCA_all.html")
    fig11.write_html(file = output_path +"/fig11a_PCA_all.html", auto_open = False)
    if fig_show: 
        fig11_loadings.show()
    with open(output_path +"/fig11b_Loadings_all.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11_loadings))
    #pyo.plot(fig11_loadings, filename = output_path +"/fig11b_Loadings_all.html")
    fig11_loadings.write_html(file = output_path +"/fig11b_Loadings_all.html", auto_open = False)

################################################################################################
    # Fig 12 PCA on raw data (only plotted if we have more than one raw file)

    
    if nr_rawfiles > 1:
        feature_list = ["RT_duration", 
        "total_num_ms1",
        "total_num_ms2", 
        "RT_TIC_Q_000-025",
        "RT_TIC_Q_025-050",
        "RT_TIC_Q_050-075",
        "RT_TIC_Q_075-100",
        "RT_MS1_Q_000-025",
        "RT_MS1_Q_025-050",
        "RT_MS1_Q_050-075",
        "RT_MS1_Q_075-100",
        "RT_MS2_Q_000-025",
        "RT_MS2_Q_025-050", 
        "RT_MS2_Q_050-075",
        "RT_MS2_Q_075-100", 
        "MS1-TIC-Change-Q2",
        "MS1-TIC-Change-Q3", 
        "MS1-TIC-Change-Q4",
        "MS1-TIC-Q2",
        "MS1-TIC-Q3",
        "MS1-TIC-Q4", 
        "MS1_Freq_Max",
        "MS1_Density_Q1",
        "MS1_Density_Q2", 
        "MS1_Density_Q3",
        "MS2_Freq_Max",
        "MS2_Density_Q1",
        "MS2_Density_Q2", 
        "MS2_Density_Q3",
        "MS2_PrecZ_1", 
        "MS2_PrecZ_2", 
        "MS2_PrecZ_3",
        "MS2_PrecZ_4",
        "MS2_PrecZ_5",
        "MS2_PrecZ_more", 
        "accumulated-MS1_TIC", 
        "accumulated-MS2_TIC"]


        ### for now, just check if the items in featurelist are really in the data frame and if not, skip them
        valid_features = [feature for feature in feature_list if feature in df.columns]
        
        df_pl12 = df[valid_features]

        df_pl12_norm = pd.DataFrame(StandardScaler().fit_transform(df_pl12)) 

        #perform PCA
        pca = PCA(n_components=2)

        principalComponents = pca.fit_transform(df_pl12_norm)

        col = range(1,(principalComponents.shape[1]+1))
        col = ['pca'+ str(y) for y in col]
        principalDf = pd.DataFrame(data = principalComponents, columns = col)

        principalDf["t_scaled"] = t_scaled
        principalDf["raw_file"] = df["filename"]

        # Explained variance in axis labels
        pca_var = pca.explained_variance_ratio_
        pca_var1 = round(pca_var[0]*100, 1)
        pca_var2 = round(pca_var[1]*100, 1)

        label_x = "PC1 (" + str(pca_var1) + "%)"
        label_y = "PC2 (" + str(pca_var2) + "%)"

        if use_group:
            principalDf["group"] = group
            fig12 = px.scatter(principalDf, x = "pca1", y = "pca2", color = "group",
                        color_continuous_scale="bluered", title = "PCA on raw data",
                        hover_name="raw_file", hover_data=["pca1", "pca2"],
                            labels={
                            "pca1": label_x,
                            "pca2": label_y,
                            "group": "Group"
                        })
        else:
            principalDf["t_scaled"] = t_scaled
            fig12 = px.scatter(principalDf, x = "pca1", y = "pca2", color = "t_scaled",
                        color_continuous_scale="bluered", title = "PCA on raw data",
                        hover_name="raw_file", hover_data=["pca1", "pca2"],
                            labels={
                            "pca1": label_x,
                            "pca2": label_y,
                            "t_scaled": "timestamp"
                        })
            
        # Table and plot with feature loadings (weights of the variables in the PCA)
        loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_features)
        loadings["length"] = np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2)
        loadings["variable"] = loadings.index
        fig12_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (raw data)", 
            hover_name="variable", hover_data=["PC1", "PC2"],)
    else:
        fig12 = go.Figure()
        fig12.add_annotation(
            x=0.5,
            y=0.5,
            text="PCA cannot be computed using only one sample!",
            showarrow=False,
            font=dict(size=14)
        )
        fig12.update_layout(
            width=600,
            height=400,
            title="Empty Plot"
        )
        fig12_loadings = fig12
    if fig_show:
        fig12.show()
    with open(output_path +"/fig12a_PCA_raw.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig12)) 
    #pyo.plot(fig12, filename = output_path +"/fig12_PCA_raw.html")
    fig12.write_html(file = output_path +"/fig12_PCA_raw.html", auto_open = False)
    if fig_show: 
        fig12_loadings.show()
    with open(output_path +"/fig12b_Loadings_raw.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig12_loadings))
    #pyo.plot(fig12_loadings, filename = output_path +"/fig12b_Loadings_raw.html")
    fig12_loadings.write_html(file = output_path +"/fig12b_Loadings_raw.html", auto_open = False)
        

#################################################################################################
    ### Fig 13: Ion Maps (one for each raw file)
    ionmap_df = []
    for index in df.index:
        tmp = dict(RT = ast.literal_eval(df["ms2_rt_array"].iloc[index]),
            MZ = ast.literal_eval(df["ms2_mz_array"].iloc[index]),
            TIC = ast.literal_eval(df["ms2_tic_array"].iloc[index]),
            filename = [df["filename"].iloc[index]]*len(ast.literal_eval(df["ms2_rt_array"].iloc[index])))
        tmp = pd.DataFrame(tmp)
        ionmap_df.append(tmp)
    ionmap_df2 = pd.concat(ionmap_df)

    # create folder to store ion maps
    if not os.path.exists(output_path + "/fig13_ionmaps"):
        os.makedirs(output_path + "/fig13_ionmaps")

    for file in df["filename"]:
        ionmap_df2_tmp = ionmap_df2[ionmap_df2["filename"] == file]
        fig13_tmp = px.density_contour(ionmap_df2_tmp, x="RT", y="MZ", title = file, nbinsx = 50, nbinsy = 50)
        fig13_tmp.update_traces(contours_coloring="fill", contours_showlabels = True)
        #if fig_show: 
            #fig12_tmp.show()
        with open(output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig13_tmp))
        #pyo.plot(fig13_tmp, filename = output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".html")
        fig13_tmp.write_html(file = output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".html", auto_open = False)

    if fig_show:
        fig13_tmp.show()


################################################################################################
    ### Figure 14: Pump Pressure


    if "THERMO_pump_pressure_bar_x_axis_____pickle_zlib" in df.columns:

        ### TODO: generate empty plot if this column is missing
        ### TODO: why is this plot also misisng for EXII???
        ### extract data from compressed columns and put them into long format
        x = []
        y = []
        fn = []
        for index in df.index:
            if pd.isnull(df["THERMO_pump_pressure_bar_x_axis_____pickle_zlib"].iloc[index]) \
                or pd.isnull(df["THERMO_pump_pressure_bar_y_axis_____pickle_zlib"].iloc[index]) :
                continue
            x_locally = unbase64_uncomp_unpickle(df["THERMO_pump_pressure_bar_x_axis_____pickle_zlib"].iloc[index])
            y_locally = unbase64_uncomp_unpickle(df["THERMO_pump_pressure_bar_y_axis_____pickle_zlib"].iloc[index])

            # TODO Add Bruker pump pressure!

            if x_locally is None:
                #print(f"x is None at {index} => {df['filename'].iloc[index]}")
                continue
            if y_locally is None:
                #print(f"y is None at {index} => {df['filename'].iloc[index]}")
                continue
            if len(x_locally) != len(y_locally):
                raise ValueError("x and y does not have same length")


            # With more than 10000 datapoints plotting the data
            # leads to unnecessary delay. Interpolating 10000 datapoints is usually enough.
            if len(x_locally) > 10000:
                samples = int(len(x_locally) / 10000)
                # Explictly adding the last datapoint to make sure we cover rounding errors when calculating `sample`
                x_locally = [x_locally[i] for i in range(0, len(x_locally), samples)] + x_locally[-1:]
                y_locally = [y_locally[i] for i in range(0, len(y_locally), samples)] + y_locally[-1:]
            x += x_locally
            y += y_locally
            fn += [df["filename"].iloc[index]] * len(x_locally)

        pp_df2 = pd.DataFrame({
            "filename": fn,
            "x": x,
            "y": y
        })
        
    else: 
        pp_df2 = pd.DataFrame()
        
    if not pp_df2.empty:
        fig14 = px.line(pp_df2, x="x", y="y", color = "filename", title = "Pump Pressure")
        fig14.update_traces(line=dict(width=0.5))
        fig14.update_yaxes(exponentformat="E") 
        fig14.update_layout(width = int(1000), height = int(1000), 
                            xaxis_title = "Time (min)", 
                            yaxis_title = "Pump pressure")
        
    else:     
        
        fig14 = go.Figure()
        fig14.add_annotation(
            x=0.5,
            y=0.5,
            text="No Pump Pressure data available!",
            showarrow=False,
            font=dict(size=14)
        )
        fig14.update_layout(
            width=600,
            height=400,
            title="Empty Plot"
        )
        
        
    if fig_show:
        fig14.show()
    with open(output_path +"/fig14_Pump_pressure.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig14))
    #pyo.plot(fig14, filename = output_path +"/fig14_Pump_pressure.html")
    fig14.write_html(file = output_path +"/fig14_Pump_pressure.html", auto_open = False)

################################################################################################
    # Figure 15_XXX: visualize all THERMO_LOG and THERMO_EXTRA (without a filter)


    for column_header in df.columns:
        if column_header.startswith("THERMO_LOG_") or column_header.startswith("THERMO_EXTRA_"):
            ### extract data from compressed columns and put them into long format
            x = [] # x-axis Scan_StartTime_zlib
            y = [] # y-axis THERMO HEADER
            fn = [] # filename

            for index in df.index:
                column_display_header = column_header

                if pd.isnull(df[column_header].iloc[index]):
                    # Skip, there is no info available
                    continue

                y_locally = unbase64_uncomp_unpickle(df[column_header].iloc[index])
                x_locally = unbase64_uncomp_unpickle(df["THERMO_Scan_StartTime_____pickle_zlib"].iloc[index])  # All of THERMO_EXTRA and THERMO_LOG are defined over the Retention tims / Scan StartTime
                
                # Keep only values for MS1 spectra  (e.g. for Lock Mass Correction or Ion Injection Time)
                if any(x in column_header for x in ["Ion Injection Time", "LM Correction", "LM m/z-Correction"]):
                    mslevel = unbase64_uncomp_unpickle(df["THERMO_Scan_msLevel_____pickle_zlib"].iloc[index])
                    x_locally = [x for x,y in zip(x_locally, mslevel) if y == 1]
                    y_locally = [float(x) for x,y in zip(y_locally, mslevel) if y == 1]
                    column_display_header = column_header.split("_____")[0] + " (MS1 Level filtered)_____" + column_header.split("_____")[1]

                x += [float(_x) for _x in x_locally]
                y += [float(_y) for _y in y_locally]
                fn += [df["filename"].iloc[index]] * len(x_locally)


            local_df = pd.DataFrame({
                "filename": fn,
                "x": x,
                "y": y
            })

            column_title = column_display_header.split("_____")[0]
            if not local_df.empty:
                fig15 = px.line(local_df, x="x", y="y", color = "filename", title = column_title)
                fig15.update_traces(line=dict(width=0.5))
                fig15.update_yaxes(exponentformat="E") 
                fig15.update_layout(width = int(1000), height = int(1000), 
                                    xaxis_title = "Time (min)", 
                                    yaxis_title = column_title)

            else: 
                fig15 = go.Figure()
                fig15.add_annotation(
                    x=0.5,
                    y=0.5,
                    text="No '{}' available!".format(column_title),
                    showarrow=False,
                    font=dict(size=14)
                )
                fig15.update_layout(
                    width=600,
                    height=400,
                    title="Empty Plot"
                )

            os.makedirs(output_path + os.sep + "THERMO_PLOTS_FIG15", exist_ok=True)
            if fig_show:
                fig15.show()
            with open(output_path + os.sep + "THERMO_PLOTS_FIG15" + os.sep + "{}.json".format(re.sub('\W+','', column_title)), "w") as json_file:
                json_file.write(plotly.io.to_json(fig15))
            fig15.write_html(file = output_path + os.sep + "THERMO_PLOTS_FIG15" + os.sep + "{}.html".format(re.sub('\W+','', column_title)), auto_open = False)


# %%
