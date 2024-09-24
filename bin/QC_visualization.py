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
#import ast
import os
#import base64
#import pickle
#import zlib
from datetime import datetime, timezone
#import math
#import pathlib

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
import h5py

pio.renderers.default = "png"


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hdf5_folder", help="Folder containing hdf5 files, one for each processed raw file", default=None)
    parser.add_argument("-group", help="List of the experimental group (comma-separated).", default=None)
    parser.add_argument("-output", help="Output folder for the plots as json files.", default = "graphics")
    parser.add_argument("-tic_overlay_offset", help = "Offset for TIC overlay plots", default = 0)
    parser.add_argument("-fig_show", help = "Show figures, e.g. for debugging?", default = False, action = "store_true")
    parser.add_argument("-spikeins", help = "Whether to analyse spike-ins", default = False, action = "store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

####################################################################################################
    # parameters

    ### get hdf5 files
    hdf5_folder = args.hdf5_folder
    hdf5_files = [f for f in os.listdir(hdf5_folder) if f.endswith('.hdf5')] 
    hdf5_file_names = [re.sub(r'\.hdf5$', '', f) for f in hdf5_files] # file names without file ending
    nr_rawfiles = len(hdf5_files)

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
    analyse_spikeins = args.spikeins


############################################################################################
# Feature list: features for initial table 0
    feature_list = ["filename", # will be filled in later as name of the hdf5 file
                    "timestamp", # will be converted later to human-readable format
                    "number_ungrouped_proteins",
                    "number_proteins",
                    "number_filtered_peptides",
                    "number_filtered_psms",
                    "total_num_ms1",
                    "total_num_ms2",
                    "Total_Ion_Current_Max",
                    "Base_Peak_Intensity_Max",
                    "Total_Ion_Current_Max_Up_To_105",
                    "Base_Peak_Intensity_Max_Up_To_105", 
                    "total_num_ident_features",
                    "total_num_features",
                    "RT_duration", 
                    "accumulated-MS1_TIC",
                    "accumulated-MS2_TIC",
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
                    "MS2_PrecZ_Unknown",
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
                    # "THERMO_TUNE_ '''Ion Transfer Tube Temperature (+ or +-)'''"
                    # "THERMO_TUNE_ '''Ion Transfer Tube Temperature (-)'''"
                    # "THERMO_TUNE_ '''Vaporizer Temp. (+ or +-)'''"
                    # "THERMO_TUNE_ '''Vaporizer Temp. (-)'''"
                ]



####################################################################################################
    # build pandas data frame with all columns for table 0
    # for each hdf5 file fill a new row with the corresponding data
    df = pd.DataFrame(columns=list(feature_list))
    for file in hdf5_files:
        hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
        new_row = pd.Series(dtype='float64')
        for feature in feature_list:
            if feature in hdf5_tmp:
                new_row[feature] = hdf5_tmp[feature][:]
                if len(new_row[feature]) == 1:
                    new_row[feature] = new_row[feature][0]
            else:
                new_row[feature] = None
        df = pd.concat([df, pd.DataFrame([new_row], columns=new_row.index)]).reset_index(drop=True)

    ## add file names
    df["filename"] = hdf5_file_names
    ## convert timestamp to something human-readable
    x = [datetime.fromtimestamp(x, timezone.utc) for x in df["timestamp"]]
    df["timestamp"] = x

    
    ### get info for each file if it is a Bruker or Thermo file
    is_bruker = []
    is_thermo = []
    add_thermo_headers = []
    add_bruker_headers = []

    for file in hdf5_files:
        hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
        keys_tmp = list(hdf5_tmp.keys())
        
        if any(keys.startswith("THERMO") for keys in hdf5_tmp.keys()):
            is_thermo.extend([True])
            is_bruker.extend([False])
            add_thermo_headers.extend([key for key in keys_tmp if key.startswith("THERMO")])
        if any(keys.startswith("BRUKER") for keys in hdf5_tmp.keys()):
            is_thermo.extend([False])
            is_bruker.extend([True])
            add_bruker_headers.extend([key for key in keys_tmp if key.startswith("BRUKER")])

    ### remove duplicates
    add_thermo_headers = set(add_thermo_headers)
    add_bruker_headers = set(add_bruker_headers)

    ############################################################################################
    # If spike-ins are analyzed, add them to the table
    if analyse_spikeins:
        spike_columns = [key for key in hdf5_tmp.keys() if re.match("SPIKE", key)]
        df_spikes = pd.DataFrame(columns=list(spike_columns))
        for file in hdf5_files:
            hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
            new_row = pd.Series(dtype='float64')
            for feature in spike_columns:
                if feature in hdf5_tmp:
                    new_row[feature] = hdf5_tmp[feature][:]
                    if len(new_row[feature]) == 1:
                        new_row[feature] = new_row[feature][0]
                else:
                    new_row[feature] = None
            df_spikes = pd.concat([df_spikes, pd.DataFrame([new_row], columns=new_row.index)]).reset_index(drop=True)
        df_table0 = pd.concat([df.loc[:,"filename":"Base_Peak_Intensity_Max_Up_To_105"], df_spikes, df.loc[:,"total_num_ident_features":]], axis=1).reindex(df.index)
    else: 
        df_table0 = df.copy()
        
    df_table0.to_csv(output_path + "/00_table_summary.csv", index = False)        

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
    fig3.write_html(file = output_path +"/fig3_barplot_features.html", auto_open = False)


####################################################################################################
    ## Figure 4: TIC Overlay as Lineplot
    
    tic_df = []
    i = 0
    for file in hdf5_files:
        hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
        tic = dict(TIC = hdf5_tmp["ms1_tic_array"][:],
                RT = hdf5_tmp["ms1_rt_array"][:],
                filename=[hdf5_file_names[i]]*len(hdf5_tmp["ms1_tic_array"][:]))
        tic = pd.DataFrame(tic)
        tic_df.append(tic)
        i += 1
    tic_df2 = pd.concat(tic_df)
        
    fig4 = px.line(tic_df2, x="RT", y="TIC", color = "filename", title = "TIC overlay")
    fig4.update_traces(line=dict(width=0.5))
    fig4.update_yaxes(exponentformat="E") 
    if fig_show:
        fig4.show()
    with open(output_path +"/fig4_TIC_overlay.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig4))
    fig4.write_html(file = output_path +"/fig4_TIC_overlay.html", auto_open = False)
    
    ### TODO: test the offset functionality
    
    '''
    offset_ratio = float(args.tic_overlay_offset)

    tic_df = []
    for index in df.index:
        tic = dict(TIC=df["ms1_tic_array"].iloc[index],
            RT=df["ms1_rt_array"].iloc[index],
            filename=[df["filename"].iloc[index]]*len(df["ms1_tic_array"].iloc[index]))
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
    '''

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
    fig8.write_html(file = output_path +"/fig8_barplot__precursor_chargestate.html", auto_open = False)

################################################################################################
    # Figure 9: PSM charge states (of identified spectra)
    
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
            width=1500,
            height=1000,
            title="Empty Plot"
        )
    if fig_show:
        fig9.show()
    with open(output_path +"/fig9_barplot_PSM_chargestate.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig9))
    fig9.write_html(file = output_path +"/fig9_barplot_PSM_chargestate.html", auto_open = False)

################################################################################################
    # Figure 10: Missed cleavages of PSMs
    
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
            width=1500,
            height=1000,
            title="Empty Plot"
        )  
    
    if fig_show:
        fig10.show()
    with open(output_path +"/fig10_barplot_PSM_missedcleavages.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig10))
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
            
        fig11.update_layout(width = int(1500), height = int(1000))
            
        # Table and plot with feature loadings (weights of the variables in the PCA)
        loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_features)
        loadings["length"] = np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2)
        loadings["variable"] = loadings.index
        fig11_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (all data)", 
            hover_name="variable", hover_data=["PC1", "PC2"],)
        fig11_loadings.update_layout(width = int(1500), height = int(1000))
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
            width=1500,
            height=1000,
            title="Empty Plot"
        )
        fig11_loadings = fig11
    if fig_show:
        fig11.show()
    with open(output_path +"/fig11a_PCA_all.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11))
    fig11.write_html(file = output_path +"/fig11a_PCA_all.html", auto_open = False)
    if fig_show: 
        fig11_loadings.show()
    with open(output_path +"/fig11b_Loadings_all.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11_loadings))
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
        fig12.update_layout(width = int(1500), height = int(1000))
            
            
        # Table and plot with feature loadings (weights of the variables in the PCA)
        loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=valid_features)
        loadings["length"] = np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2)
        loadings["variable"] = loadings.index
        fig12_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (raw data)", 
            hover_name="variable", hover_data=["PC1", "PC2"],)
        fig12_loadings.update_layout(width = int(1500), height = int(1000))
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
            width=1500,
            height=1000,
            title="Empty Plot"
        )
        fig12_loadings = fig12
    if fig_show:
        fig12.show()
    with open(output_path +"/fig12a_PCA_raw.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig12)) 
    fig12.write_html(file = output_path +"/fig12_PCA_raw.html", auto_open = False)
    if fig_show: 
        fig12_loadings.show()
    with open(output_path +"/fig12b_Loadings_raw.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig12_loadings))
    fig12_loadings.write_html(file = output_path +"/fig12b_Loadings_raw.html", auto_open = False)
        

#################################################################################################
    ### Fig 13: Ion Maps (one for each raw file)

    ionmap_df = []
    i = 0
    for file in hdf5_files:
        hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
        tmp = dict(RT = hdf5_tmp["ms2_rt_array"][:],
                MZ = hdf5_tmp["ms2_mz_array"][:],
                filename=[hdf5_file_names[i]]*len(hdf5_tmp["ms2_rt_array"][:]))
        tmp = pd.DataFrame(tmp)
        ionmap_df.append(tmp)
        i += 1
    ionmap_df2 = pd.concat(ionmap_df)

    # create folder to store ion maps
    if not os.path.exists(output_path + "/fig13_ionmaps"):
        os.makedirs(output_path + "/fig13_ionmaps")

    for file in df["filename"]:
        ionmap_df2_tmp = ionmap_df2[ionmap_df2["filename"] == file]
        fig13_tmp = px.density_contour(ionmap_df2_tmp, x="RT", y="MZ", title = file, nbinsx = 50, nbinsy = 50)
        fig13_tmp.update_traces(contours_coloring="fill", contours_showlabels = True)
        fig13_tmp.update_layout(width = 1500, height = 1000, 
                    xaxis_title = "Retention Time (seconds)", 
                    yaxis_title = "m/z")
        with open(output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig13_tmp))
        fig13_tmp.write_html(file = output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".html", auto_open = False)

    if fig_show:
        fig13_tmp.show()


################################################################################################
    ### Figure 14: Pump Pressure

    ### only for Thermo, Bruker pump pressure is plotted as one of the extra plots in figure 15
    if any(is_thermo):
        pump_df = []
        i = 0
        for file in hdf5_files:
            #file = hdf5_files[i]
            hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
            if "THERMO_pump_pressure_bar_x_axis" in hdf5_tmp.keys() and "THERMO_pump_pressure_bar_y_axis" in hdf5_tmp.keys():
                x = hdf5_tmp["THERMO_pump_pressure_bar_x_axis"][:]
                y = hdf5_tmp["THERMO_pump_pressure_bar_y_axis"][:]
                
                # use not more than 10000 data points. If data has more than 10000 data points, take every nth data point        
                if len(x) > 10000: 
                    samples = int(len(x) / 10000)
                    x = [x[i] for i in range(0, len(x), samples)] 
                    y = [y[i] for i in range(0, len(y), samples)] 
                    
                pump_df_tmp = dict(filename = hdf5_file_names[i],
                                    x = x,
                                    y = y)
                pump_df_tmp = pd.DataFrame(pump_df_tmp)
                
                pump_df.append(pump_df_tmp)
            i += 1
            

        pump_df2 = pd.DataFrame()
        if type(x) == list:
            pump_df2 = pd.concat(pump_df)
        elif type(x) == pd.DataFrame:  # if only one raw file has pump pressure data available
            pump_df2 = pump_df

        ### TODO: what if none of the files has pump pressure data available?

        fig14 = px.line(pump_df2, x="x", y="y", color = "filename", title = "Pump Pressure")
        fig14.update_traces(line=dict(width=0.5))
        fig14.update_yaxes(exponentformat="E") 
        fig14.update_layout(width = 1500, height = 1000, 
                            xaxis_title = "Time (min)", 
                            yaxis_title = "Pump pressure")
        fig14.write_html(file = output_path +"/fig14_Pump_pressure.html", auto_open = False)
        
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
            width=1500,
            height=1000,
            title="Empty Plot"
        )
        
        
    if fig_show:
        fig14.show()
    with open(output_path +"/fig14_Pump_pressure.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig14))
    fig14.write_html(file = output_path +"/fig14_Pump_pressure.html", auto_open = False)




################################################################################################
    # Figure 15_XXX: visualize all THERMO_LOG and THERMO_EXTRA data

    os.makedirs(output_path + os.sep + "THERMO_PLOTS_FIG15", exist_ok=True) 
    if any(is_thermo):
        for header in add_thermo_headers:
            if header.startswith("THERMO_LOG_") or header.startswith("THERMO_EXTRA_"):  # this exludes the pump pressure data, MS level and scan start time
                ### extract data from compressed columns and put them into long format
                x = [] # x-axis Scan_StartTime_zlib
                y = [] # y-axis THERMO HEADER
                fn = [] # filename
                
                i = 0
                for file in hdf5_files:
                #for index in df.index:
                    hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')
                    #column_display_header = header

                    if header not in hdf5_tmp.keys():
                        # Skip, there is no info available for this hdf5 file
                        i += 1
                        continue
                    
                    #if type(df[header].iloc[index]) is float \
                    #    and math.isnan(df[header].iloc[index]):
                    #    continue

                    y_tmp = hdf5_tmp[header][:]
                    x_tmp = hdf5_tmp["THERMO_Scan_StartTime"][:]  # All of THERMO_EXTRA and THERMO_LOG are defined over the Retention tims / Scan StartTime
                    
                    # Keep only values for MS1 spectra  (e.g. for Lock Mass Correction or Ion Injection Time)
                    display_header = header
                    if any(x in header for x in ["Ion Injection Time", "LM Correction", "LM m/z-Correction"]):
                        mslevel = hdf5_tmp["THERMO_Scan_msLevel"][:]
                        x_tmp = [x for x,y in zip(x_tmp, mslevel) if y == 1]
                        y_tmp = [float(x) for x,y in zip(y_tmp, mslevel) if y == 1]
                        display_header = header + " (MS1 Level filtered)" ## add info about MS1 level filtering to the plot title

                    x += [float(_x) for _x in x_tmp]
                    y += [float(_y) for _y in y_tmp]
                    fn += [hdf5_file_names[i]] * len(x_tmp)
                    i += 1


                df_tmp = pd.DataFrame({
                    "filename": fn,
                    "x": x,
                    "y": y
                })
                
                #print(df_tmp.head())
                #column_title = column_display_header.split("_____")[0]
                if not df_tmp.empty:
                    fig15 = px.line(df_tmp, x = "x", y = "y", color = "filename", title = display_header)
                    fig15.update_traces(line = dict(width = 0.5))
                    fig15.update_yaxes(exponentformat="E") 
                    fig15.update_layout(width = int(1500), height = int(1000), 
                                        xaxis_title = "Time (min)", 
                                        yaxis_title = display_header)
                else: 
                    fig15 = go.Figure()
                    fig15.add_annotation(
                        x=0.5,
                        y=0.5,
                        text="No '{}' available!".format(display_header),
                        showarrow=False,
                        font=dict(size=14)
                    )
                    fig15.update_layout(
                        width=1500,
                        height=1000,
                        title="Empty Plot"
                    )

                #os.makedirs(output_path + os.sep + "THERMO_PLOTS_FIG15", exist_ok=True)
                if fig_show:
                    fig15.show()
                with open(output_path + os.sep + "THERMO_PLOTS_FIG15" + os.sep + "{}.json".format(re.sub('\W+','', display_header)), "w") as json_file:
                    json_file.write(plotly.io.to_json(fig15))
                fig15.write_html(file = output_path + os.sep + "THERMO_PLOTS_FIG15" + os.sep + "{}.html".format(re.sub('\W+','', display_header)), auto_open = False)
        
        
        
        
        
    ################################################################################################
        # Figure 16_XXX: visualize all additional BRUKER data

    os.makedirs(output_path + os.sep + "BRUKER_PLOTS_FIG16", exist_ok=True)
    if any(is_bruker): 
        for header in add_bruker_headers:
           
            if header == "BRUKER_Time": # Time, will be needed as x-axis in all plots
                continue
            if header == "BRUKER_MsMsType":  # MsMsType codes for DIA/DDA for example. Doesn't need to be plotted.
                continue
            if header.startswith("BRUKER_pump_pressure_bar"):  # Skip Pump Pressure
                continue
            
            ### extract data from compressed columns and put them into long format
            x = [] # x-axis Scan_StartTime_zlib
            y = [] # y-axis BRUKER HEADER
            fn = [] # filename
            
            i = 0
            for file in hdf5_files:
                hdf5_tmp = h5py.File(hdf5_folder + "/" + file,'r')

                if header not in hdf5_tmp.keys():
                    # Skip, there is no info available for this hdf5 file
                    i += 1
                    continue
                
                y_tmp = hdf5_tmp[header][:]
                x_tmp = hdf5_tmp["BRUKER_Time"][:]  # All BRUKER variables are defined over the time
                
                display_header = header

                x += [float(_x) for _x in x_tmp]
                y += [float(_y) for _y in y_tmp]
                fn += [hdf5_file_names[i]] * len(x_tmp)
                i += 1


            df_tmp = pd.DataFrame({
                "filename": fn,
                "x": x,
                "y": y
            })
        
            if not df_tmp.empty:
                fig16 = px.line(df_tmp, x="x", y="y", color = "filename", title = display_header)
                fig16.update_traces(line=dict(width=0.5))
                fig16.update_yaxes(exponentformat="E") 
                fig16.update_layout(width = int(1500), height = int(1000), 
                                    xaxis_title = "Time (min)", 
                                    yaxis_title = display_header)
            else: 
                fig16 = go.Figure()
                fig16.add_annotation(
                    x=0.5,
                    y=0.5,
                    text="No '{}' available!".format(display_header),
                    showarrow=False,
                    font=dict(size=14)
                )
                fig16.update_layout(
                    width=1500,
                    height=1000,
                    title="Empty Plot"
                )

            if fig_show:
                fig16.show()
            with open(output_path + os.sep + "BRUKER_PLOTS_FIG16" + os.sep + "{}.json".format(re.sub('\W+','', display_header)), "w") as json_file:
                json_file.write(plotly.io.to_json(fig16))
            fig16.write_html(file = output_path + os.sep + "BRUKER_PLOTS_FIG16" + os.sep + "{}.html".format(re.sub('\W+','', display_header)), auto_open = False)

        
      

# %%
