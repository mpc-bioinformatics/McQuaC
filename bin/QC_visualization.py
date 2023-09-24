#!/bin/env python

# %%

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
import plotly.offline as pyo
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-csv_file", help="CSV file with QC data.")
    parser.add_argument("-group", help="List of the experimental group (comma-separated).", default=None)
    parser.add_argument("-output", help="Output folder for the plots as json files.", default = "graphics")
    parser.add_argument("-tic_overlay_offset", help = "Offset for TIC overlay plots", default = 0)
    parser.add_argument("-fig_show", help = "Show figures, e.g. for debugging?", default = False, action = "store_true")
    parser.add_argument("-isa", help = "Is this ISA QC?", default = False, action = "store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()
    print(args)

### TODO: plots auch als plotly-html ausgeben

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

####################################################################################################
    # short function to extract infos from compressed columns
    def unbase64_uncomp_unpickle(x: bytes) -> list:
        unb64 = base64.b64decode(x)
        uncomp = zlib.decompress(unb64)
        return pickle.loads(uncomp)


    if isa:
       ############################################################################################
        # Table 0: table with overviw over QC measures. 
        feature_list = ["filename",
                        "timestamp",
                        "number_ungrouped_proteins",
                        "number_proteins", 
                        "number_filtered_peptides", 
                        "number_filtered_pms", # TODO: typo
                        "total_num_ms1", 
                        "total_num_ms2",                      
                        "Total_Ion_Current_Max",
                        "Base_Peak_Intensity_Max",
                        "Total_Ion_Current_Max_Up_To_105",
                        "FIXED_MPCSPIKE1_PEP_old-isa01_MZ_457.2834_RT_2439",
                        "FIXED_MPCSPIKE2_PEP_old-isa02_MZ_895.9493_RT_4212",
                        "FIXED_MPCSPIKE3_PEP_old-isa03_MZ_966.088_RT_5925",
                        "FIXED_MPCSPIKE4_PEP_GEPAAAAAPEAGASPVEK_MZ_815.9118_RT_1890",
                        "FIXED_MPCSPIKE5_PEP_NLVVGDETTSSLR_MZ_700.8664_RT_2611",
                        "FIXED_MPCSPIKE6_PEP_LQPGDIGIYR_MZ_571.3156_RT_3067",
                        "FIXED_MPCSPIKE7_PEP_VVVLPSGALQISR_MZ_674.913_RT_3995",
                        "FIXED_MPCSPIKE8_PEP_YPGAYYIFQIK_MZ_685.8654_RT_4595",
                        "FIXED_MPCSPIKE9_PEP_NIPTVNENLENYYLEVNQLEK_MZ_848.7618_RT_5490",
                        "IDENT_MPCSPIKE1_COUNT",
                        "IDENT_MPCSPIKE1_DELTA_RT",
                        "IDENT_MPCSPIKE1_PEP_old-isa01_MZ_457.2834_RT_DELTA",
                        "IDENT_MPCSPIKE2_COUNT",
                        "IDENT_MPCSPIKE2_DELTA_RT",
                        "IDENT_MPCSPIKE2_PEP_old-isa02_MZ_895.9493_RT_DELTA",
                        "IDENT_MPCSPIKE3_COUNT",
                        "IDENT_MPCSPIKE3_DELTA_RT",
                        "IDENT_MPCSPIKE3_PEP_old-isa03_MZ_966.088_RT_DELTA",
                        "IDENT_MPCSPIKE4_COUNT",
                        "IDENT_MPCSPIKE4_DELTA_RT",
                        "IDENT_MPCSPIKE4_PEP_GEPAAAAAPEAGASPVEK_MZ_815.9118_RT_DELTA",
                        "IDENT_MPCSPIKE5_COUNT",
                        "IDENT_MPCSPIKE5_DELTA_RT",
                        "IDENT_MPCSPIKE5_PEP_NLVVGDETTSSLR_MZ_700.8664_RT_DELTA",
                        "IDENT_MPCSPIKE6_COUNT",
                        "IDENT_MPCSPIKE6_DELTA_RT",
                        "IDENT_MPCSPIKE6_PEP_LQPGDIGIYR_MZ_571.3156_RT_DELTA",
                        "IDENT_MPCSPIKE7_COUNT",
                        "IDENT_MPCSPIKE7_DELTA_RT",
                        "IDENT_MPCSPIKE7_PEP_VVVLPSGALQISR_MZ_674.913_RT_DELTA",
                        "IDENT_MPCSPIKE8_COUNT",
                        "IDENT_MPCSPIKE8_DELTA_RT",
                        "IDENT_MPCSPIKE8_PEP_YPGAYYIFQIK_MZ_685.8654_RT_DELTA",
                        "IDENT_MPCSPIKE9_COUNT",
                        "IDENT_MPCSPIKE9_DELTA_RT",
                        "IDENT_MPCSPIKE9_PEP_NIPTVNENLENYYLEVNQLEK_MZ_848.7618_RT_DELTA"
                    ]

        x = [datetime.datetime.utcfromtimestamp(x) for x in df["timestamp"]] # convert timestamp to datetime
        df_isatable0 = df[feature_list]
        df_isatable0 = df_isatable0.loc[:,:].copy()
        df_isatable0.loc[:,"timestamp"] = x
        df_isatable0.to_csv(output_path + "/isatable0_summary.csv")

        
        #   Figure 1: Number of protein groups and accessions
        df_pl1_isa = df[["filename", "number_proteins", "number_ungrouped_proteins"]]
        df_pl1_long = df_pl1_isa.melt(id_vars = ["filename"])

        fig1_isa = px.bar(df_pl1_long, x="filename", y="value", color="variable", barmode = "group", 
                    title = "Number of protein groups and accessions")
        fig1_isa.update_yaxes(exponentformat="none") 
        if fig_show: 
            fig1_isa.show()

        with open(output_path +"/isafig1_barplot_proteingroups.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig1_isa))
        pyo.plot(fig1_isa, filename = output_path +"/isafig1_barplot_proteingroups.html")


        # Figure 2: Number of peptides
        df_pl2_isa = df[["filename", "number_filtered_peptides"]]

        fig2_isa = px.bar(df_pl2_isa, x="filename", y="number_filtered_peptides",
                    title = "Number of peptides")
        fig2_isa.update_yaxes(exponentformat="none") 

        with open(output_path +"/isafig2_barplot_peptides.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig2_isa))
        pyo.plot(fig2_isa, filename = output_path +"/isafig2_barplot_peptides.html")

        
        # Figure 3: Number of PSMs
        #### TODO: has to be tested, when PIA output is finished!
        df_pl3_isa = df[["filename", "number_filtered_pms"]] # TODO: typo in psms

        fig3_isa = px.bar(df_pl3_isa, x="filename", y="number_filtered_pms", 
                    title = "Number of PSMs")
        fig3_isa.update_yaxes(exponentformat="none") 

        with open(output_path +"/isafig3_barplot_PSMs.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig3_isa))
        pyo.plot(fig3_isa, filename = output_path +"/isafig3_barplot_PSMs.html")
        
            
    else:  # normal QC

        ############################################################################################
        # Table 0: table with overviw over QC measures. 
        feature_list = ["filename",
                        "timestamp",
                        "total_num_ms1",
                        "total_num_ms2",
                        "number_filtered_pms", # TODO: typo
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
                        "Z1",  # TODO: changed column name
                        "Z2", 
                        "Z3", 
                        "Z4", 
                        "Z5", 
                        "missed_0",
                        "missed_1",
                        "missed_2",
                        "missed_3",
                        "missed_more"
                    ]
        
        x = [datetime.datetime.utcfromtimestamp(x) for x in df["timestamp"]] # convert timestamp to datetime
        df_table0 = df.loc[:,:].copy()
        df_table0.loc[:,"timestamp"] = x
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
        pyo.plot(fig1, filename = output_path +"/fig1_barplot_MS1_MS2.html")
        


    ################################################################################################
        # Figure 2: Barplot for number of PSMs, peptides, proteins
        df_pl2 = df[["filename", "number_filtered_psms", "number_filtered_peptides", "number_proteins", "number_ungrouped_proteins"]]
        df_pl2_long = df_pl2.melt(id_vars = ["filename"])
        fig2 = px.bar(df_pl2_long, x="filename", y="value", color="variable", barmode = "group", 
                    title = "Number of filtered PSMs, filtered peptides, filtered protein groups and accessions")
        fig2.update_yaxes(exponentformat="none") 
        if fig_show: 
            fig2.show()
        with open(output_path +"/fig2_barplot_PSMs_peptides_proteins.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig2))
        pyo.plot(fig2, filename = output_path +"/fig2_barplot_PSMs_peptides_proteins.html")

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
        pyo.plot(fig3, filename = output_path +"/fig3_barplot_features.html")


    ####################################################################################################
        ## Figure 4: TIC Overlay as Lineplot
        offset_ratio = float(args.tic_overlay_offset)

        #MS1_TICs = df[["ms1_tic_array", "ms1_rt_array"]]
        tic_df = []
        for index in df.index:
            tic = dict(TIC=ast.literal_eval(df["ms1_tic_array"].iloc[index]),
                RT=ast.literal_eval(df["ms1_rt_array"].iloc[index]),
                filename=[df["filename"].iloc[index]]*len(ast.literal_eval(df["ms1_tic_array"].iloc[index])))
            tic = pd.DataFrame(tic)
            tic_df.append(tic)
        tic_df2 = pd.concat(tic_df)
        #print(tic_df2)

        #TIC_list2 = pd.concat(TIC_list)
        offset = max(tic_df2["TIC"])
        #print(offset)

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
        #fig3.update_layout(width = int(1000), height = int(1000))
        if fig_show:
            fig4.show()
        with open(output_path +"/fig3_TIC_overlay.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig4))
        pyo.plot(fig4, filename = output_path +"/fig4_TIC_overlay.html")

   
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
        pyo.plot(fig5, filename = output_path +"/fig5_barplot_TIC_quartiles.html")

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
        pyo.plot(fig6, filename = output_path +"/fig6_barplot_MS1_TIC_quartiles.html")   

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
        pyo.plot(fig7, filename = output_path +"/fig7_barplot_MS2_TIC_quartiles.html")

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
        pyo.plot(fig8, filename = output_path +"/fig8_barplot__precursor_chargestate.html")

    ################################################################################################
        # Figure 9: PSM charge states (of identified spectra)
        #### TODO: relative statt absolute Zahlen!
        df_pl9 = df[["filename", 'psm_charge1', 'psm_charge2', 'psm_charge3', 'psm_charge4', 'psm_charge5', 'psm_charge_more']]
        df_pl9_long = df_pl9.melt(id_vars = ["filename"])
        fig9 = px.bar(df_pl9_long, x="filename", y="value", color="variable", title = "Charge states of PSMs")
        fig9.update_xaxes(tickangle=-90)
        if fig_show:
            fig9.show()
        with open(output_path +"/fig9_barplot_PSM_chargestate.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig9))
        pyo.plot(fig9, filename = output_path +"/fig9_barplot_PSM_chargestate.html")

    ################################################################################################
        # Figure 10: Missed cleavages of PSMs
        ### TODO: relative statt absolute Zahlen!
        df_pl10 = df[["filename", 'psm_missed_0', 'psm_missed_1', 'psm_missed_2', 'psm_missed_3', 'psm_missed_more']]
        df_pl10_long = df_pl10.melt(id_vars = ["filename"])
        fig10 = px.bar(df_pl10_long, x="filename", y="value", color="variable", title = "Number of missed cleavages for PSMs")
        fig10.update_xaxes(tickangle=-90)
        if fig_show:
            fig10.show()
        with open(output_path +"/fig10_barplot_PSM_missedcleavages.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig10))
        pyo.plot(fig10, filename = output_path +"/fig10_barplot_PSM_missedcleavages.html")


    #################################################################################################
        #### Preparations for PCA
        # If no groups are given, the points in the PCA are coloured by the timestamp.
        # The oldest raw file is coloured blue, the newest red.
        # Depending on the timestamp, the other ones are coloured on a gradient in between these two colours.

        timestamps = df[["timestamp"]].values.flatten().tolist()
        mintime = min(timestamps)
        maxtime = max(timestamps)
        #print(mintime)
        #print(maxtime)

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

            df_pl11 = df[feature_list]
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
                            color_continuous_scale="bluered", title = "PCA on raw data",
                            hover_name="raw_file", hover_data=["pca1", "pca2"],
                                labels={
                                "pca1": label_x,
                                "pca2": label_y,
                                "t_scaled": "timestamp"
                            })
                
            # Table and plot with feature loadings (weights of the variables in the PCA)
            loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=feature_list)
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
            fig11_loadings = fig10
        if fig_show:
            fig11.show()
        with open(output_path +"/fig11a_PCA_all.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig11))
        pyo.plot(fig11, filename = output_path +"/fig11a_PCA_all.html")
        if fig_show: 
            fig11_loadings.show()
        with open(output_path +"/fig11b_Loadings_all.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig11_loadings))
        pyo.plot(fig11_loadings, filename = output_path +"/fig11b_Loadings_all.html")

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


            df_pl12 = df[feature_list]

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
            loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=feature_list)
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
        pyo.plot(fig12, filename = output_path +"/fig12_PCA_raw.html")
        if fig_show: 
            fig12_loadings.show()
        with open(output_path +"/fig12b_Loadings_raw.json", "w") as json_file:
            json_file.write(plotly.io.to_json(fig12_loadings))
        pyo.plot(fig12_loadings, filename = output_path +"/fig12b_Loadings_raw.html")
            

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
            fig13_tmp = px.density_contour(ionmap_df2_tmp, x="RT", y="MZ", title = file)
            fig13_tmp.update_traces(contours_coloring="fill", contours_showlabels = True)
            #if fig_show: 
                #fig12_tmp.show()
            with open(output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".json", "w") as json_file:
                json_file.write(plotly.io.to_json(fig13_tmp))
            pyo.plot(fig13_tmp, filename = output_path +"/fig13_ionmaps/fig13_ionmap_" + file + ".html")

        if fig_show:
            fig13_tmp.show()


    ################################################################################################
        ### Figure 14: Pump Pressure

        ### extract data from compressed columns and put them into long format
        x = []
        y = []
        fn = []
        for index in df.index:
            if pd.isnull(df["THERMO_pump_preasure_bar_x_axis"].iloc[index]):
                # Skip, there is no Pump pressure available
                continue
            x_locally = unbase64_uncomp_unpickle(df["THERMO_pump_preasure_bar_x_axis"].iloc[index])
            y_locally = unbase64_uncomp_unpickle(df["THERMO_pump_preasure_bar_y_axis"].iloc[index])

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

        if not pp_df2.empty:
            fig14 = px.line(pp_df2, x="x", y="y", color = "filename", title = "Pump Pressure")
            fig14.update_traces(line=dict(width=0.5))
            fig14.update_yaxes(exponentformat="E") 
            fig14.update_layout(width = int(1000), height = int(1000), 
                                xaxis_title = "Time (min)", 
                                yaxis_title = "Pump pressure")
            if fig_show:
                fig14.show()
            with open(output_path +"/fig14_Pump_pressure.json", "w") as json_file:
                json_file.write(plotly.io.to_json(fig14))
            pyo.plot(fig14, filename = output_path +"/fig14_Pump_pressure.html")
                    

    ################################################################################################
        # Figure 15: Ion Injection time

        ### extract data from compressed columns and put them into long format
        x = [] # x-axis Scan_StartTime_zlib
        y = [] # y-axis Ion_Injection_Time_pickle_zlib
        fn = [] # filename
        
        for index in df.index:

            if pd.isnull(df["THERMO_Ion_Injection_Time_pickle_zlib"].iloc[index]):
                # Skip, there is no Ion Injection Time available
                continue

            y_locally = unbase64_uncomp_unpickle(df["THERMO_Ion_Injection_Time_pickle_zlib"].iloc[index])
            y_locally = [y_locally[i] for i in range(0, len(y_locally), 2)]  ### TODO: Nur zum Testen weil komischweise doppelte Werte vorhanden
            x_locally = unbase64_uncomp_unpickle(df["THERMO_Scan_StartTime_zlib"].iloc[index])
            mslevel = unbase64_uncomp_unpickle(df["THERMO_Scan_msLevel_zlib"].iloc[index])
              
            # keep only values for MS1 spectra
            x_locally = [x for x,y in zip(x_locally, mslevel) if y == 1]
            y_locally = [float(x) for x,y in zip(y_locally, mslevel) if y == 1]
  
            x += x_locally
            y += y_locally
            fn += [df["filename"].iloc[index]] * len(x_locally)
            
            
        ionInjTime_df = pd.DataFrame({
            "filename": fn,
            "x": x,
            "y": y
        })

        if not ionInjTime_df.empty:
            fig15 = px.line(ionInjTime_df, x="x", y="y", color = "filename", title = "Ion Injection Time")
            fig15.update_traces(line=dict(width=0.5))
            fig15.update_yaxes(exponentformat="E") 
            fig15.update_layout(width = int(1000), height = int(1000), 
                                xaxis_title = "Time (min)", 
                                yaxis_title = "Ion Injection Time (ms)")
            if fig_show:
                fig15.show()
            with open(output_path +"/fig15_Ion_Injection_Time.json", "w") as json_file:
                json_file.write(plotly.io.to_json(fig15))
            pyo.plot(fig15, filename = output_path +"/fig15_Ion_Injection_Time.html")


    ################################################################################################
        # Figure 16: Lock Mass Correction

        x = [] # x-axis Scan_StartTime_zlib
        y = [] # y-axis Ion_Injection_Time_pickle_zlib
        fn = [] # filename

        for index in df.index:

            if pd.isnull(df["THERMO_LM_m_z_Correction_pickle_zlib"].iloc[index]):
                continue

            y_locally = unbase64_uncomp_unpickle(df["THERMO_LM_m_z_Correction_pickle_zlib"].iloc[index])
            x_locally = unbase64_uncomp_unpickle(df["THERMO_Scan_StartTime_zlib"].iloc[index])
            y_locally = [y_locally[i] for i in range(0, len(y_locally), 2)]  ### TODO: Nur zum Testen weil komischweise doppelte Werte vorhanden
            
            # With more than 10000 datapoints plotting the data
            # leads to unnecessary delay. Interpolating 10000 datapoints is usually enough.
            if len(x_locally) > 10000:
                samples = int(len(x_locally) / 10000)
                # Explictly adding the last datapoint to make sure we cover rounding errors when calculating `sample`
                x_locally = [x_locally[i] for i in range(0, len(x_locally), samples)] + x_locally[-1:]
                y_locally = [y_locally[i] for i in range(0, len(y_locally), samples)] + y_locally[-1:]
            
            
            y_locally = [float(x) for x in y_locally]
            x += x_locally
            y += y_locally
            fn += [df["filename"].iloc[index]] * len(x_locally)
            
        LMCorr_df = pd.DataFrame({
        "filename": fn,
        "x": x,
        "y": y
        })

        if not LMCorr_df.empty:
            fig16 = px.line(LMCorr_df, x="x", y="y", color = "filename", title = "Lock Mass Correction")
            fig16.update_traces(line=dict(width=0.5))
            fig16.update_yaxes(exponentformat="E") 
            fig16.update_layout(width = int(1000), height = int(1000),
                                xaxis_title = "Time (min)", 
                                yaxis_title = "LM m/z-Correction (ppm)")
            if fig_show:
                fig16.show()
            with open(output_path +"/fig16_Lock_Mass_Correction.json", "w") as json_file:
                json_file.write(plotly.io.to_json(fig16))
            pyo.plot(fig16, filename = output_path +"/fig16_Lock_Mass_Correction.html")
            
            
            

