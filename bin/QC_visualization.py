#!/bin/env python

# %%

# std imports
import argparse
from typing import Optional
import math
from io import BytesIO
from pathlib import Path
import zipfile
import ast

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
    parser.add_argument("-csv_file", help="CSV file with QC data.")
    parser.add_argument("-group", help="List of the experimental group (comma-separated).", default=None)
    parser.add_argument("-output", help="Output folder for the plots as json files.", default = "graphics")
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    #############
    # parameters

    ### use HCC dataset for testing
    csv_file = args.csv_file#.split(",")
    #["OEI06439", "OEI06441","OEI06443","OEI06445","OEI06447","OEI06449","OEI06453","OEI06455","OEI06459","OEI06461","OEI06463","OEI06465","OEI06467",
    #    "OEI06469","OEI06471","OEI06473","OEI06475","OEI06477","OEI06481","OEI06483","OEI06485","OEI06487","OEI06489","OEI06491","OEI06493","OEI06495","OEI06497",
    #    "OEI06499","OEI06504","OEI06506","OEI06508","OEI06510","OEI06518","OEI06520","OEI06522","OEI06526","OEI06528","OEI06530"]

    ### grouping
    #x = np.array([["HCC"],["C"]])
    ### If use_groups = False, PCA plots are coloured by timestamp. If True, PCA plots are coloured by group.
    if (args.group is None):
        use_group = False
    else:
        use_group = True
        group = np.array(args.group.split(",")) #np.repeat(x, 19)
    
    ### folder to save the plots as json files
    output_path = args.output    # "graphics/"

    df = pd.read_csv(csv_file)

    # Figure 1: Barplot for total number of MS1 and MS2 spectra
    df_pl1 = df[["filename", "total_num_ms1", "total_num_ms2"]]
    df_pl1_long = df_pl1.melt(id_vars = ["filename"])
    fig1 = px.bar(df_pl1_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Total number of MS1 and MS2 spectra")
    fig1.update_yaxes(exponentformat="none") 
    fig1.show()
    with open(output_path +"/fig1_barplot_MS1_MS2.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig1))


    # Figure 2: Barplot for number of PSMs, peptides, proteins and identified features
    ### TODO: spaeter eine Grafik mit filtered_PSMs, filtered_peptides und filtered_proteingroups und eine getrennte mit
    ###       nr_features und identified_nr_features
    ### TODO: testen, wenn PIA Teil fertig.
    '''
    df_pl2a = df[["filename", "number-filtered-psms", "number-filtered-peptides", "number-filtered-protein-groups"]]
    df_pl2a_long = df_pl2a.melt(id_vars = ["filename"])
    fig2a = px.bar(df_pl2a_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of filtered PSMs, filtered peptides and filtered protein groups")
    fig2a.update_yaxes(exponentformat="E") 
    fig2a.show()
    with open(output_path +"/fig2a_barplot_PSMs_peptides_proteins.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig2a))
    '''
    df_pl2b = df[["filename", "total_num_features", "total_num_ident_features"]]
    df_pl2b_long = df_pl2b.melt(id_vars = ["filename"])
    fig2b = px.bar(df_pl2b_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of features and identified features")
    fig2b.update_yaxes(exponentformat="none") 
    fig2b.show()
    with open(output_path +"/fig2b_barplot_features.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig2b))



    ## Figure 3: TIC Overlay as Lineplot
    MS1_TICs = df[["ms1_tic_array", "ms1_rt_array"]]
    tic_df = []
    for index in df.index:
        tic = dict(TIC=ast.literal_eval(df["ms1_tic_array"].iloc[index]),
            RT=ast.literal_eval(df["ms1_rt_array"].iloc[index]),
            filename=[df["filename"].iloc[index]]*len(ast.literal_eval(df["ms1_tic_array"].iloc[index])))
        tic = pd.DataFrame(tic)
        tic_df.append(tic)
    tic_df = pd.concat(tic_df)

    fig3 = px.line(tic_df, x="RT", y="TIC", color = "filename", title = "TIC overlay")
    fig3.update_yaxes(exponentformat="E") 
    fig3.show()
    with open(output_path +"/fig3_TIC_overlay.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig3))


    # Figure 4: Barplot TIC quartiles
    df_pl4 = df[["filename", 'RT_TIC_Q_000-025', 'RT_TIC_Q_025_050', 'RT_TIC_Q_050-075', 'RT_TIC_Q_075-100']]
    df_pl4_long = df_pl4.melt(id_vars = ["filename"])
    fig4 = px.bar(df_pl4_long, x = "filename", y = "value", color = "variable", title = "Quartiles of TIC over retention time")
    fig4.show()
    with open(output_path +"/fig4_barplot_TIC_quartiles.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig4))


    # Figure 5: Barplot MS1 TIC quartiles
    df_pl5 = df[["filename", 'RT_MS1_Q_000-025', 'RT_MS1_Q_025_050', 'RT_MS1_Q_050-075', 'RT_MS1_Q_075-100']]
    df_pl5_long = df_pl5.melt(id_vars = ["filename"])
    fig5 = px.bar(df_pl5_long, x="filename", y="value", color="variable", title = "Quartiles of MS1 over retention time")
    fig5.show()
    with open(output_path +"/fig5_barplot_MS1_TIC_quartiles.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig5))


    # Figure 6: Barplot MS2 TIC quartiles
    df_pl6 = df[["filename", 'RT_MS2_Q_000-025', 'RT_MS2_Q_025_050', 'RT_MS2_Q_050-075', 'RT_MS2_Q_075-100']]
    df_pl6_long = df_pl6.melt(id_vars = ["filename"])
    fig6 = px.bar(df_pl6_long, x="filename", y="value", color="variable", title = "Quartiles of MS2 over retention time")
    fig6.show()
    with open(output_path +"/fig6_barplot_MS2_TIC_quartiles.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig6))

    # Figure 7: Precursor charge states
    df_pl7 = df[["filename", 'MS2_PrecZ1', 'MS2_PrecZ2', 'MS2_PrecZ3', 'MS2_PrecZ4', 'MS2_PrecZ5', 'MS2_PrecZ_more']]
    df_pl7_long = df_pl7.melt(id_vars = ["filename"])
    fig7 = px.bar(df_pl7_long, x="filename", y="value", color="variable", title = "Charge states of precursors")
    fig7.show()
    with open(output_path +"/fig7_barplot_precursor_chargestate.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig7))


    # Figure 8: PSM charge states (of identified spectra)
    #### TODO: muss getestet werden, wenn PIA Teil fertig ist!
    #df_pl8 = df[["filename", 'psmZ1', 'psmZ2', 'psmZ3', 'psmZ4', 'psmZ5']]
    #df_pl8_long = df_pl8.melt(id_vars = ["filename"])
    #fig8 = px.bar(df_pl8_long, x="filename", y="value", color="variable", title = "Charge states of PSMs")
    #with open(output_path +"/fig8_barplot_PSM_chargestate.json", "w") as json_file:
    #    json_file.write(plotly.io.to_json(fig8))

    # Figure 9: Missed cleavages of PSMs
    ### TODO: muss getestet werden, wenn PIA Teil fertig ist
    #df_pl9 = df[["filename", 'psm-missed-0', 'psm-missed-1', 'psm-missed-2', 'psm-missed-3', 'psm-missed-more']]
    #df_pl9_long = df_pl9.melt(id_vars = ["filename"])
    #fig9 = px.bar(df_pl9_long, x="filename", y="value", color="variable", title = "Number of missed cleavages for PSMs")
    #fig9.show()
    #with open(output_path +"/fig9_barplot_missedcleavages_PSMs.json", "w") as json_file:
    #    json_file.write(plotly.io.to_json(fig9))

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
        t_scaled_tmp = (t - mintime)/(maxtime-mintime)*100 
        t_scaled.append(t_scaled_tmp)


    # Fig 10 PCA on all data
    #### TODO: Testen, wenn PIA Teil fertig ist
    '''
    feature_list = ["RT_duration", 
    "total_num_ms1",
    "total_num_ms2", 
    "RT_TIC_Q_000-025",
    "RT_TIC_Q_025_050",
    "RT_TIC_Q_050-075",
    "RT_TIC_Q_075-100",
    "RT_MS1_Q_000-025",
    "RT_MS1_Q_025_050",
    "RT_MS1_Q_050-075",
    "RT_MS1_Q_075-100",
    "RT_MS2_Q_000-025",
    "RT_MS2_Q_025_050", 
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
    "MS2_PrecZ1", 
    "MS2_PrecZ2", 
    "MS2_PrecZ3",
    "MS2_PrecZ4",
    "MS2_PrecZ5",
    "MS2_PrecZ_more", 
    "accumulated-MS1_TIC", 
    "accumulated-MS2_TIC",
    "total_num_ident_features",
    "num_feature_charge_1",
    "num_feature_charge_2",
    "num_feature_charge_3",
    "num_feature_charge_4",
    "num_feature_charge_5", 
    "psmZ1",
    "psmZ2", 
    "psmZ3", 
    "psmZ4", 
    "psmZ5"]

    df_pl10 = df[feature_list]

    ### Scale data (substract mean and divide by standard deviation) before calculating the PCA
    df_pl10_norm = pd.DataFrame(StandardScaler().fit_transform(df_pl10)) 

    # calculate PCA (2 components)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(df_pl10_norm)

    col = range(1,(principalComponents.shape[1]+1))
    col = ['pca'+ str(y) for y in col]
    principalDf = pd.DataFrame(data = principalComponents, columns = col)
    principalDf["raw_file"] = join_df["filename"]

    # Explained variance in axis labels
    pca_var = pca.explained_variance_ratio_
    pca_var1 = round(pca_var[0]*100, 1)
    pca_var2 = round(pca_var[1]*100, 1)
    label_x = "PC1 (" + str(pca_var1) + "%)"
    label_y = "PC2 (" + str(pca_var2) + "%)"

    # Generate the PCA plot as a Scatterplot
    if use_group:
        principalDf["group"] = group
        fig10 = px.scatter(principalDf, x = "pca1", y = "pca2", color = "group",
                    color_continuous_scale="bluered", title = "PCA on all data",
                    hover_name="raw_file", hover_data=["pca1", "pca2"],
                        labels={
                        "pca1": label_x,
                        "pca2": label_y,
                        "group": "Group"
                    })
    else:
        principalDf["t_scaled"] = t_scaled
        fig10 = px.scatter(principalDf, x = "pca1", y = "pca2", color = "t_scaled",
                    color_continuous_scale="bluered", title = "PCA on all data",
                    hover_name="raw_file", hover_data=["pca1", "pca2"],
                        labels={
                        "pca1": label_x,
                        "pca2": label_y,
                        "t_scaled": "timestamp"
                    })
        
    #fig10.show()

    with open(output_path +"/fig10a_PCA_all.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig10))

    # %%
    # Table and plot with feature loadings (weights of the variables in the PCA)
    loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=feature_list)
    loadings["length"] = np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2)
    loadings["variable"] = loadings.index

    fig10_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (all data)", 
    hover_name="variable", hover_data=["PC1", "PC2"],)
    #fig10_loadings.show()

    with open(output_path +"/fig10b_Loadings_all.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig10_loadings))
    '''

    # Fig 11 PCA on raw data
    feature_list = ["RT_duration", 
    "total_num_ms1",
    "total_num_ms2", 
    "RT_TIC_Q_000-025",
    "RT_TIC_Q_025_050",
    "RT_TIC_Q_050-075",
    "RT_TIC_Q_075-100",
    "RT_MS1_Q_000-025",
    "RT_MS1_Q_025_050",
    "RT_MS1_Q_050-075",
    "RT_MS1_Q_075-100",
    "RT_MS2_Q_000-025",
    "RT_MS2_Q_025_050", 
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
    "MS2_PrecZ1", 
    "MS2_PrecZ2", 
    "MS2_PrecZ3",
    "MS2_PrecZ4",
    "MS2_PrecZ5",
    "MS2_PrecZ_more", 
    "accumulated-MS1_TIC", 
    "accumulated-MS2_TIC"]


    df_pl11 = df[feature_list]

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
                    color_continuous_scale="bluered", title = "PCA on raw data",
                    hover_name="raw_file", hover_data=["pca1", "pca2"],
                        labels={
                        "pca1": label_x,
                        "pca2": label_y,
                        "t_scaled": "timestamp"
                    })
        
    fig11.show()
    with open(output_path +"/fig11a_PCA_raw.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11))


    # Table and plot with feature loadings (weights of the variables in the PCA)
    loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=feature_list)
    loadings["length"] = np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2)
    loadings["variable"] = loadings.index
    fig11_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (raw data)", 
    hover_name="variable", hover_data=["PC1", "PC2"],)
    fig11_loadings.show()
    with open(output_path +"/fig11b_Loadings_raw.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11_loadings))


