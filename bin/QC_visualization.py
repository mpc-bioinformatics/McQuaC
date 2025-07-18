#!/usr/bin/env python

# %%
from collections import defaultdict
from pathlib import Path
import re
# std imports
import argparse
import os
from datetime import datetime, timezone
from typing import Dict, Tuple, List, Any

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
import matplotlib.pyplot as plt

pio.renderers.default = "png"


def get_dataset_types(hdf: h5py.File) -> Tuple[Tuple[str], Tuple[str], Tuple[str]]:
    """
    Get the types of the datasets in the hdf5 file

    Parameters
    ----------
    hdf : h5py.File
        The hdf5 file

    Returns
    -------
    Tuple[Tuple[str], Tuple[str], Tuple[str]]
        A tuple of three tuples. The first tuple contains the names of the single value datasets, 
        the second tuple contains the names of the array datasets 
        and the third tuple contains the names of the datafram datasets.
    """
    single_value_ids = []
    array_value_ids = []
    dataframe_ids = []

    for key in hdf.keys():
        if isinstance(hdf[key], h5py.Dataset):
            if hdf[key].shape[0] == 1:
                ### TODO!: bug if first file has only a single PSM, then psm-error is put in the false category
                if key == "LOCAL:05|filtered_psms_ppm_error":
                    continue
                single_value_ids.append(key)
            else:
                array_value_ids.append(key)
        elif isinstance(hdf[key], h5py.Group):
            dataframe_ids.append(key)

    return (tuple(single_value_ids), tuple(array_value_ids), tuple(dataframe_ids))

def get_dataframe_of_single_values(
        hdfs: List[h5py.File], single_value_ids: Tuple[str]
) -> pd.DataFrame:
    """
    Get a dataframe of the single values of the hdf5 files

    Parameters
    ----------
    hdfs : List[h5py.File]
        List of hdf5 files
    single_value_ids : Tuple[str]
        Tuple of the ids of the single values

    Returns
    -------
    pd.DataFrame
        The dataframe with the single values
    """
    single_value_data: Dict[str, List[Any]] = {"filename": [Path(hdf.filename).stem for hdf in hdfs]}
    for sv_id in single_value_ids:
        short_name = sv_id.split("|")[-1]
        single_value_data[short_name] = [hdf[sv_id][0] for hdf in hdfs]
    return pd.DataFrame(single_value_data)

def get_array_values(
        hdfs: List[h5py.File], array_value_ids: Tuple[str]
) -> Dict[str, Dict[str, List[Any]]]:
    """
    Get a dictionary of the arrays of the hdf5 files

    Parameters
    ----------
    hdfs : List[h5py.File]
        List of hdf5 files
    array_value_ids : Tuple[str]
        Tuple of the ids of the array values

    Returns
    -------
    Dict[str, Dict[str, List[Any]]]
        The dictionary with the array values, e.g. {filename: {array_name: [array_values]}}
    """
    array_value_data: Dict[str, Dict[str, List[Any]]] = defaultdict(dict)
    for hdf in hdfs:
        filename = Path(hdf.filename).stem
        for array_id in array_value_ids:
            short_name = array_id.split("|")[-1]
            array_value_data[filename][short_name] = hdf[array_id][:]

    return array_value_data


def get_dataframes_values(
        hdfs: List[h5py.File], dataframe_ids: Tuple[str]
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Get a dictionary of the dataframes (matrices) of the hdf5 files

    Parameters
    ----------
    hdfs : List[h5py.File]
        List of hdf5 files
    dataframe_ids : Tuple[str]

    Returns
    -------
    Dict[str, Dict[str, List[Any]]]
        The dictionary with the dataframes, e.g. {filename: {dataframe_name: pd.DataFrame}}
    """

    dataframes: Dict[str, Dict[str, List[Any]]] = defaultdict(dict)
    for hdf in hdfs:
        filename = Path(hdf.filename).stem
        for df_id in dataframe_ids:#
            short_name = df_id.split("|")[-1]
            column_order = hdf[df_id].attrs["column_order"].split("|")
            df = pd.DataFrame(dict(hdf[df_id]))
            df = df[column_order]
            dataframes[filename][short_name] = df

    return dataframes


def assemble_result_table(
    metric_list: List[str],
    hdf5_file_names: List[str],
    single_values: pd.DataFrame,
    single_value_ids_short: List[str],
    array_values: Dict[str, Dict[str, List[Any]]],
    array_value_ids_short: List[str],
    dataframes: Dict[str, Dict[str, pd.DataFrame]],
    dataframe_ids_short: List[str],
    spikein_columns: List[str] = ["Maximum_Intensity", "RT_at_Maximum_Intensity", "PSMs", "Delta_to_expected_RT"], 
    RT_unit: str = "sec"
    ) -> pd.DataFrame:
    """
    Assemble the result table
    
    Parameters
    ----------
    
    metric_list : List[str]
        List of metrics to be included in the table
    hdf5_file_names : List[str]
        List of hdf5 file names
    single_values : pd.DataFrame
        DataFrame with the single values
    single_value_ids_short : List[str]
        List of the short names of the single values
    array_values : Dict[str, Dict[str, List[Any]]]
        Dictionary with the array values
    array_value_ids_short : List[str]
        List of the short names of the array values
    dataframes : Dict[str, Dict[str, pd.DataFrame]]
        Dictionary with the dataframes
    dataframe_ids_short : List[str]
        List of the short names of the dataframes
    spikein_columns : List[str]
        List of the columns of the spike-in dataframes that should end up in the result table.
        Must be a subset of ["Maximum_Intensity", "RT_at_Maximum_Intensity", "PSMs", "Delta_to_expected_RT"].
    RT_unit : str
        Unit of the retention time, either sec for seconds or min for minutes.
        
    Returns
    -------
    pd.DataFrame
        The result table
    """

    df_table = pd.DataFrame()
    for metric in metric_list:
        if metric == "filename":
            df_table["filename"] = hdf5_file_names
       
        elif metric in single_value_ids_short:
            if metric == "timestamp":
                ## convert timestamp to something human-readable
                x = [datetime.fromtimestamp(x, timezone.utc) for x in single_values["timestamp"]]
                df_table[metric] = x
                df_table['timestamp'] = df_table['timestamp'].dt.tz_localize(None)
            else:
                df_table[metric] = single_values[metric]
        
        ### match base peak intensity and total ion current with upper time limit
        elif metric in ["base_peak_intensity_max_up_to_", "total_ion_current_max_up_to_"]:
                matching_metrics = [element for element in single_value_ids_short if metric in element]
                df_table[matching_metrics] = single_values[matching_metrics]
        
        elif metric in array_value_ids_short:
            df_tmp = pd.DataFrame()
            for file in hdf5_file_names:
                values = array_values[file][metric]
                if metric == "RT_range":
                    columns = [metric + "_" + suffix for suffix in ["min", "max"]]
                else: 
                    columns = [metric + "_" + str(i) for i in range(1, len(values)+1)]
                df_tmp_tmp = pd.DataFrame(columns = columns)
                df_tmp_tmp.loc[0] = values
                df_tmp = pd.concat([df_tmp, df_tmp_tmp], axis = 0)
                
      
                
            df_tmp.reset_index(drop=True, inplace=True)
            ### transform time to minutes if necessary
            if metric == "RT_range":
                if RT_unit == "min":
                    df_tmp["RT_range_min"] = df_tmp["RT_range_min"]/60
                    df_tmp["RT_range_max"] = df_tmp["RT_range_max"]/60
            
            df_table = pd.concat([df_table, df_tmp], axis = 1)
        
        elif metric in dataframe_ids_short:
            if metric == "spike_in_metrics":
                df_table_spike = pd.DataFrame()
                for file in hdf5_file_names:
                    spike_data = dataframes[file][metric]
                    
                    if RT_unit == "min":
                        spike_data["RT_at_Maximum_Intensity"] = spike_data["RT_at_Maximum_Intensity"]/60
                        spike_data["Delta_to_expected_RT"] = spike_data["Delta_to_expected_RT"]/60
                    
                    spike_in_list = spike_data['Spike-in'].astype(str).tolist()
                    #spike_name = [s.split("_")[1] for s in spike_in_list]
                    mz = [s.split("_")[3] for s in spike_in_list]
                    spike_data["mz"] = ["MZ_" + x for x in mz]
                    column_names_spike = ["SPIKE_" + x for x in spike_data["mz"].astype(str).tolist()]
                    
                    ## for each spike-in, extract the data
                    df_tmp = pd.DataFrame()  ## data frame for each spike-in and each file
                    for index, row in spike_data.iterrows():
                        column_names = [column_names_spike[index] + "_" + x for x in spikein_columns]    
                        df_tmp_tmp = pd.DataFrame(columns = column_names)
                        df_tmp_tmp.loc[0] = spike_data.loc[index, spikein_columns].values

                        df_tmp = pd.concat([df_tmp, df_tmp_tmp], axis = 1) ## add columns
                        
                    df_table_spike = pd.concat([df_table_spike, df_tmp], axis = 0)
                    df_table_spike.reset_index(drop=True, inplace=True)
                    
                df_table = pd.concat([df_table, df_table_spike], axis = 1)

                    
            else:
                df_tmp = pd.DataFrame()
                for file in hdf5_file_names:
                    df_tmp_tmp = dataframes[file][metric]
                    df_tmp = pd.concat([df_tmp, df_tmp_tmp], axis = 0)
                columns = [metric + "_" + str(col) for col in df_tmp.columns]
                df_tmp.reset_index(drop=True, inplace=True)
                df_tmp.columns = columns
                df_table = pd.concat([df_table, df_tmp], axis = 1)
     
        else:  
            df_table[metric] = np.nan ### fill column with NaNs if metric is not in the hdf5 files
            
    return df_table






def check_if_file_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    if os.path.isfile(s):
        return s
    else:
        raise Exception("File '{}' does not exists".format(s))


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hdf5_files", type=check_if_file_exists, nargs="+", help = "hdf5 files which are used for visualization as string separated by whitespace", default = None)
    parser.add_argument("-output", help="Output folder for the plots as json files.", default = "graphics")
    parser.add_argument("-output_table_type", help="Type of output table (one of csv, tsv or xlsx)", default = "csv")
    parser.add_argument("-spikeins", help = "Whether to analyse spike-ins", default = False, action = "store_true")
    parser.add_argument("-group", help="List of the experimental group (comma-separated).", default=None)  ### TODO: input table with group information
    parser.add_argument("-RT_unit", help="Unit of the retention time, either sec for seconds or min for minutes.", default = "sec")
    parser.add_argument("-fig_show", help = "Show figures, e.g. for debugging?", default = False, action = "store_true")
    parser.add_argument("-output_column_order", help = "Order of columns in the output table", default = "", type = str)
    parser.add_argument("-spikein_columns", help = "Columns of the spike-in dataframes that should end up in the result table", default = "Maximum_Intensity,RT_at_Maximum_Intensity,PSMs,Delta_to_expected_RT", type = str)
    parser.add_argument("-height_barplots", help = "Height of the barplots in pixels", default = 700, type = int) # in pixels
    parser.add_argument("-width_barplots", help = "Width of the barplots in pixels", default = 0, type = int) # default 0: flexible width, in pixels
    parser.add_argument("-height_pca", help = "Height of the PCA plots in pixels", default = 1000, type = int) # in pixels
    parser.add_argument("-width_pca", help = "Width of the PCA plots in pixels", default = 1000, type = int) # in pixels
    parser.add_argument("-height_ionmaps", help = "Height of the ionmaps in inches", default = 10, type = int)
    parser.add_argument("-width_ionmaps", help = "Width of the ionmaps in inches", default = 10, type = int)
    return parser.parse_args()


##########################################################################################################################################################
### main function to visualize QC results

if __name__ == "__main__":
    args = argparse_setup()
    print(args)

    args.hdf5_files = sorted(args.hdf5_files) # sorts the file names alphabetically (assumes that they are all in the same folder)

    ### read in the hdf5 files
    hdf5s = [h5py.File(f, "r") for f in args.hdf5_files]

    (single_value_ids, array_value_ids, dataframe_ids) =  get_dataset_types(hdf5s[0])
    
    single_values = get_dataframe_of_single_values(hdf5s, single_value_ids)
    #single_values = single_values.sort_values("filename", ascending=True)  # Sort values by filename
    single_value_ids_short = [s.split("|")[-1] for s in single_value_ids]

    array_values = get_array_values(hdf5s, array_value_ids)
    array_value_ids_short = [s.split("|")[-1] for s in array_value_ids]

    ### remove Thermo headers with "Extracted_Log_Headers" from dataframe_ids (cannot be plotted because only single values)
    dataframe_ids = [x for x in dataframe_ids if x != "THERMO_LOG|Extracted_Log_Headers"]

    dataframes = get_dataframes_values(hdf5s, dataframe_ids)
    dataframe_ids_short = [s.split("|")[-1] for s in dataframe_ids]

            
            

####################################################################################################
    # parameters

    hdf5_file_names = single_values["filename"].values
    nr_rawfiles = len(hdf5_file_names)
    
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
    
    
    ##########################################################################################
    ### order of columns for output table
    
    if args.output_column_order == "":  ### take default column order
        metric_list = [
            "filename",
            "timestamp",
            "RT_range",
            "nr_MS1",
            "nr_MS2",
            "accumulated_MS1_TIC",
            "accumulated_MS2_TIC",
            "base_peak_intensity_max",
            "total_ion_current_max",
            "base_peak_intensity_max_up_to_",
            "total_ion_current_max_up_to_",
            "MS2_prec_charge_fraction",
            "RT_MS1_quartiles",
            "RT_MS2_quartiles",
            "RT_TIC_quartiles",
            "MS1_freq_max",
            "MS2_freq_max",
            "MS1_density_quartiles",
            "MS2_density_quartiles",
            "MS1_TIC_change_quartiles",
            "MS1_TIC_quartiles",
            "nr_PSMs",
            "nr_peptides",
            "nr_protein_groups",
            "nr_accessions",
            "PSM_charge_fractions",
            "PSM_missed_cleavage_counts",
            "nr_features",
            "nr_ident_features",
            "features_charge",
            "ident_features_charge",
            "spike_in_metrics"
        ]
    else:  ### take user-defined column order
        metric_list = args.output_column_order.split(",")
        
    df_table0 = assemble_result_table(
        metric_list = metric_list, 
        hdf5_file_names = hdf5_file_names,
        single_values = single_values,
        single_value_ids_short = single_value_ids_short,
        array_values = array_values,
        array_value_ids_short = array_value_ids_short,
        dataframes = dataframes,
        dataframe_ids_short = dataframe_ids_short,
        RT_unit = args.RT_unit, 
        spikein_columns = args.spikein_columns.split(",")
    )

    # Sort values by filename
    #df_table0 = df_table0.sort_values(by = "filename", ascending=True)  
    single_values = single_values.sort_values(by = "filename", ascending=True)


    if args.output_table_type == "csv":
        df_table0.to_csv(output_path + os.sep + "00_table_summary.csv", index = False)
    if args.output_table_type == "tsv":   
        df_table0.to_csv(output_path + os.sep + "00_table_summary.csv", index = False, sep = "\t")
    if args.output_table_type == "xlsx":
        df_table0.to_excel(output_path + os.sep + "00_table_summary.xlsx", index = False)   
    

################################################################################################
    # Figure 01: Barplot for total number of MS1 and MS2 spectra
    df_pl01 = single_values[["filename", "nr_MS1", "nr_MS2"]]
    df_pl01_long = df_pl01.melt(id_vars = ["filename"])
    fig01 = px.bar(df_pl01_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of MS1 and MS2 spectra")
    fig01.update_yaxes(exponentformat="none") 
    fig01.update_xaxes(tickangle=-90)
    fig01.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig01.update_layout(width = int(args.width_barplots))
    
    if fig_show: 
        fig01.show()
    with open(output_path + os.sep + "fig01_barplot_MS1_MS2.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig01))
    fig01.write_html(file = output_path + os.sep + "fig01_barplot_MS1_MS2.html", auto_open = False)



################################################################################################
    # Figure 02: Barplot for number of PSMs, peptides, proteins
    
    df_pl02 = single_values[["filename", "nr_PSMs", "nr_peptides", "nr_protein_groups", "nr_accessions"]]
    df_pl02_long = df_pl02.melt(id_vars = ["filename"])
    fig02 = px.bar(df_pl02_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of filtered PSMs, filtered peptides, filtered protein groups and accessions")
    fig02.update_yaxes(exponentformat="none") 
    fig02.update_xaxes(tickangle=-90)
    fig02.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig02.update_layout(width = int(args.width_barplots))
    if fig_show: 
        fig02.show()
    with open(output_path + os.sep + "fig02_barplot_PSMs_peptides_proteins.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig02))
    fig02.write_html(file = output_path + os.sep + "fig02_barplot_PSMs_peptides_proteins.html", auto_open = False)
    
    
################################################################################################
    # Figure 03: Barplot for features and identified features
    
    df_pl03 = single_values[["filename", "nr_features", "nr_ident_features"]]
    df_pl03_long = df_pl03.melt(id_vars = ["filename"])
    fig03 = px.bar(df_pl03_long, x="filename", y="value", color="variable", barmode = "group", 
                title = "Number of features and identified features")
    fig03.update_yaxes(exponentformat="none") 
    fig03.update_xaxes(tickangle=-90)
    fig03.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig03.update_layout(width = int(args.width_barplots))
    if fig_show: 
        fig03.show()
    with open(output_path + os.sep + "fig03_barplot_features.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig03))
    fig03.write_html(file = output_path + os.sep + "fig03_barplot_features.html", auto_open = False)


####################################################################################################
    ## Figure 04: TIC Overlay as Lineplot
    
    tic_df = []
    for file in hdf5_file_names:
        tic_tmp = dataframes[file]["MS1_TIC"]
        tic_tmp["filename"] = [file]*len(tic_tmp)
        tic_df.append(tic_tmp)
    tic_df2 = pd.concat(tic_df)
    
    if args.RT_unit == "min":
        tic_df2["retention_time"] = tic_df2["retention_time"]/60
    
    #tic_df2 = tic_df2.sort_values(by = ["filename", "retention_time"], ascending=True)  
    
    fig04 = px.line(tic_df2, x="retention_time", y="TIC", color = "filename", title = "TIC overlay")
    fig04.update_traces(line=dict(width=0.5))
    fig04.update_yaxes(exponentformat="E") 
    fig04.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig04.update_layout(width = int(args.width_barplots))
    if args.RT_unit == "sec":
        fig04.update_layout(xaxis_title = "Retention Time (sec)")
    elif args.RT_unit == "min":
        fig04.update_layout(xaxis_title = "Retention Time (min)")

    if fig_show:
        fig04.show()
    with open(output_path + os.sep + "fig04_MS1_TIC_overlay.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig04))
    fig04.write_html(file = output_path + os.sep + "fig04_MS1_TIC_overlay.html", auto_open = False)
    

 
#################################################################################################
    # Figure 05: Barplot TIC quartiles
    
    RT_TIC_Q_df_list = []
    for file in hdf5_file_names:
        df_tmp = pd.DataFrame()
        df_tmp["filename"] = [file]*4
        df_tmp["variable"] = ["RT_TIC_Q_" + str(i) for i in range(1,5)]
        df_tmp["value"] = array_values[file]["RT_TIC_quartiles"]
        RT_TIC_Q_df_list.append(df_tmp)
    df_pl05_long = pd.concat(RT_TIC_Q_df_list)
    #df_pl05_long = df_pl05_long.sort_values(by = "filename", ascending=True)  
   
    fig05 = px.bar(df_pl05_long, x = "filename", y = "value", color = "variable", title = "Quartiles of TIC over retention time")
    fig05.update_xaxes(tickangle=-90)
    fig05.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig05.update_layout(width = int(args.width_barplots))
    if fig_show:
        fig05.show()
    with open(output_path + os.sep + "fig05_barplot_TIC_quartiles.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig05))
    fig05.write_html(file = output_path + os.sep + "fig05_barplot_TIC_quartiles.html", auto_open = False)




################################################################################################
    # Figure 06: Barplot MS1 TIC quartiles
    
    RT_MS1_Q_df_list = []
    for file in hdf5_file_names:
        df_tmp = pd.DataFrame()
        df_tmp["filename"] = [file]*4
        df_tmp["variable"] = ["RT_MS1_Q_" + str(i) for i in range(1,5)]
        df_tmp["value"] = array_values[file]["RT_MS1_quartiles"]
        RT_MS1_Q_df_list.append(df_tmp)
    df_pl06_long = pd.concat(RT_MS1_Q_df_list)
    #df_pl06_long = df_pl06_long.sort_values(by = "filename", ascending=True)  
    
    fig06 = px.bar(df_pl06_long, x="filename", y="value", color="variable", title = "Quartiles of MS1 over retention time")
    fig06.update_xaxes(tickangle=-90)
    fig06.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig06.update_layout(width = int(args.width_barplots))
    if fig_show: 
        fig06.show()
    with open(output_path + os.sep + "fig06_barplot_MS1_TIC_quartiles.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig06))
    fig06.write_html(file = output_path + os.sep + "fig06_barplot_MS1_TIC_quartiles.html", auto_open = False)
    

################################################################################################
    # Figure 07: Barplot MS2 TIC quartiles
    
    RT_MS2_Q_df_list = []
    for file in hdf5_file_names:
        df_tmp = pd.DataFrame()
        df_tmp["filename"] = [file]*4
        df_tmp["variable"] = ["RT_MS2_Q_" + str(i) for i in range(1,5)]
        df_tmp["value"] = array_values[file]["RT_MS2_quartiles"]
        RT_MS2_Q_df_list.append(df_tmp)
    df_pl07_long = pd.concat(RT_MS2_Q_df_list)
    #df_pl07_long = df_pl07_long.sort_values(by = "filename", ascending=True)  
    
    fig07 = px.bar(df_pl07_long, x="filename", y="value", color="variable", title = "Quartiles of MS2 over retention time")
    fig07.update_xaxes(tickangle=-90)
    fig07.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig07.update_layout(width = int(args.width_barplots))
    if fig_show:
        fig07.show()
    with open(output_path + os.sep + "fig07_barplot_MS2_TIC_quartiles.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig07))
    fig07.write_html(file = output_path + os.sep + "fig07_barplot_MS2_TIC_quartiles.html", auto_open = False)
    

################################################################################################
    # Figure 08: Precursor charge states
    
    Prec_charge_df_list = []
    for file in hdf5_file_names:
        df_tmp = dataframes[file]["MS2_prec_charge_fraction"]
        df_tmp.rename(columns = {"0": "unknown", df_tmp.columns[-1]: "more"}, inplace = True)
        df_tmp_long = df_tmp.melt()
        df_tmp_long["filename"] = [file]*df_tmp.shape[1]
        Prec_charge_df_list.append(df_tmp_long)
    df_pl08_long = pd.concat(Prec_charge_df_list)
    df_pl08_long.rename(columns = {"variable": "Prec_charge", "value": "fraction"}, inplace = True)
    #df_pl08_long = df_pl08_long.sort_values(by = "filename", ascending=True)  
    
    fig08 = px.bar(df_pl08_long, x="filename", y="fraction", color="Prec_charge", title = "Charge states of precursors")
    fig08.update_xaxes(tickangle=-90)
    fig08.update_layout(height = int(args.height_barplots))
    if args.width_barplots > 0:
        fig08.update_layout(width = int(args.width_barplots))
    if fig_show:
        fig08.show()
    with open(output_path + os.sep + "fig08_barplot_precursor_charge.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig08))
    fig08.write_html(file = output_path + os.sep + "fig08_barplot_precursor_charge.html", auto_open = False)


################################################################################################
    # Figure 09: PSM charge states (of identified spectra)
    
    if ("PSM_charge_fractions" in dataframes[hdf5_file_names[0]].keys()):
        PSM_charge_df_list = []
        for file in hdf5_file_names:
            df_tmp = dataframes[file]["PSM_charge_fractions"]
            df_tmp.rename(columns = {"0": "unknown", df_tmp.columns[-1]: "more"}, inplace = True)
            df_tmp_long = df_tmp.melt()
            df_tmp_long["filename"] = [file]*df_tmp.shape[1]
            PSM_charge_df_list.append(df_tmp_long)
        df_pl09_long = pd.concat(PSM_charge_df_list)
        df_pl09_long.rename(columns = {"variable": "PSM_charge", "value": "fraction"}, inplace = True)
        #df_pl09_long = df_pl09_long.sort_values(by = "filename", ascending=True)  
        
        fig09 = px.bar(df_pl09_long, x="filename", y="fraction", color="PSM_charge", title = "Charge states of PSMs")
        fig09.update_xaxes(tickangle=-90)
        fig09.update_layout(height = int(args.height_barplots))
        if args.width_barplots > 0:
            fig09.update_layout(width = int(args.width_barplots))
    else: 
        fig09 = go.Figure()
        fig09.add_annotation(
            x=0.5,
            y=0.5,
            text="Columns are missing, no plot created!",
            showarrow=False,
            font=dict(size=14)
        )
        fig09.update_layout(
            width=1500,
            height=1000,
            title="Empty Plot"
        )
    if fig_show:
        fig09.show()
    with open(output_path + os.sep + "fig09_barplot_PSM_charge.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig09))
    fig09.write_html(file = output_path + os.sep + "fig09_barplot_PSM_charge.html", auto_open = False)


################################################################################################
    # Figure 10: Missed cleavages of PSMs (normalized between 0 and 1 in "fractions")

    if ("PSM_missed_cleavage_counts" in dataframes[hdf5_file_names[0]].keys()):
        PSM_missed_df_list = []
        for file in hdf5_file_names:
            df_tmp = dataframes[file]["PSM_missed_cleavage_counts"]
            df_tmp.rename(columns = {df_tmp.columns[-1]: "more"}, inplace = True)
            df_tmp_long = df_tmp.melt()
            df_tmp_long["filename"] = [file]*df_tmp.shape[1]
            PSM_missed_df_list.append(df_tmp_long)
        df_pl10_long = pd.concat(PSM_missed_df_list)
        df_pl10_long_perc = df_pl10_long.copy()
        df_pl10_long_perc["value"] = df_pl10_long["value"]/df_pl10_long.groupby("filename")["value"].transform("sum")
        df_pl10_long_perc.rename(columns = {"variable": "PSM_missed_cleavages", "value": "Fraction"}, inplace = True)
        #df_pl10_long_perc = df_pl10_long_perc.sort_values(by = "filename", ascending=True)  
        
        fig10 = px.bar(df_pl10_long_perc, x="filename", y="Fraction", color="PSM_missed_cleavages", title = "Fraction of missed cleavages for PSMs")
        fig10.update_xaxes(tickangle=-90)
        fig10.update_layout(height = int(args.height_barplots))
        if args.width_barplots > 0:
            fig10.update_layout(width = int(args.width_barplots))
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
    with open(output_path + os.sep + "fig10_barplot_PSM_missedcleavages.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig10))
    fig10.write_html(file = output_path + os.sep + "fig10_barplot_PSM_missedcleavages.html", auto_open = False)
    
################################################################################################
    # Fig 11 PCA on raw data (before identification)
    # (only plotted if we have more than one raw file)

    ## Calculate scaling of timestamps to colour the points in the PCA plot (percentage between min and max time):
    timestamps = single_values["timestamp"].values.flatten().tolist()
    mintime = min(timestamps)
    maxtime = max(timestamps)
    
    t_scaled = []
    for t in timestamps:
        if (mintime == maxtime):
            t_scaled_tmp = 1
        else:
            t_scaled_tmp = (t - mintime)/(maxtime-mintime)*100 
        t_scaled.append(t_scaled_tmp)

    
    if nr_rawfiles > 1:
        metric_list_PCA_raw = ["RT_range", 
        "nr_MS1",
        "nr_MS2", 
        "accumulated_MS1_TIC", 
        "accumulated_MS2_TIC",
        "base_peak_intensity_max",
        "total_ion_current_max",
        "MS2_prec_charge_fraction",
        "RT_MS1_quartiles",
        "RT_MS2_quartiles",
        "RT_TIC_quartiles",
        "MS1_freq_max",
        "MS2_freq_max",
        "MS1_density_quartiles",
        "MS2_density_quartiles",
        "MS1_TIC_change_quartiles",
        "MS1_TIC_quartiles"]

        df_pl11 = assemble_result_table(
            metric_list = metric_list_PCA_raw, 
            hdf5_file_names = hdf5_file_names,
            single_values = single_values,
            single_value_ids_short = single_value_ids_short,
            array_values = array_values,
            array_value_ids_short = array_value_ids_short,
            dataframes = dataframes,
            dataframe_ids_short = dataframe_ids_short
        )
 
        df_pl11 = df_pl11.fillna(value = 0) # impute missig values by 0
        ## scale the data before computing the PCA
        df_pl11_scaled = pd.DataFrame(StandardScaler().fit_transform(df_pl11)) 

        #perform PCA
        pca = PCA(n_components=2)

        principalComponents = pca.fit_transform(df_pl11_scaled)

        col = range(1,(principalComponents.shape[1]+1))
        col = ['pca'+ str(y) for y in col]
        principalDf = pd.DataFrame(data = principalComponents, columns = col)

        principalDf["t_scaled"] = t_scaled
        principalDf["raw_file"] = hdf5_file_names

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
        fig11.update_layout(width = int(args.width_pca), height = int(args.height_pca))
        fig11.update_traces(marker=dict(size=20))
            
            
        # Table and plot with feature loadings (weights of the variables in the PCA)
        loadings = pd.DataFrame(pca.components_.T, columns = ['PC1', 'PC2'], index = df_pl11.columns) 
        loadings.insert(0, "length", np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2))
        loadings.insert(0, "variable", loadings.index)
        loadings.sort_values("length", ascending=False, inplace=True)
        fig11_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (raw data)", 
            hover_name="variable", hover_data=["PC1", "PC2"],)
        fig11_loadings.update_layout(width = int(args.width_pca), height = int(args.height_pca))
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
        loadings = pd.DataFrame(columns = ["variable", "length", "PC1", "PC2"])
    if fig_show:
        fig11.show()
    with open(output_path + os.sep + "fig11a_PCA_raw.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11)) 
    fig11.write_html(file = output_path + os.sep + "fig11_PCA_raw.html", auto_open = False)
    if fig_show: 
        fig11_loadings.show()
    with open(output_path + os.sep + "fig11b_Loadings_raw.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11_loadings))
    fig11_loadings.write_html(file = output_path + os.sep + "fig11b_Loadings_raw.html", auto_open = False)
    
    
    ### save loadings as tables
    if args.output_table_type == "csv":
        loadings.to_csv(output_path + os.sep + "fig11c_table_loadings_raw.csv", index = False)
    if args.output_table_type == "tsv":   
        loadings.to_csv(output_path + os.sep + "fig11c_table_loadings_raw.tsv", index = False, sep = "\t")
    if args.output_table_type == "xlsx":
        loadings.to_excel(output_path + os.sep + "fig11c_table_loadings_raw.xlsx", index = False)   
       

#################################################################################################
    # Fig 12 PCA on all data (after identification)
    # (only plotted if we have more than one raw file)

    if nr_rawfiles > 1:
        metric_list_PCA_all = metric_list_PCA_raw + [
            "nr_PSMs",
            "nr_peptides",
            "nr_protein_groups",
            "nr_accessions",
            "PSM_charge_fractions",
            "PSM_missed_cleavage_counts",
            "nr_features",
            "nr_ident_features",
            "features_charge",
            "ident_features_charge"]
        
        df_pl12 = assemble_result_table(
            metric_list = metric_list_PCA_all, 
            hdf5_file_names = hdf5_file_names,
            single_values = single_values,
            single_value_ids_short = single_value_ids_short,
            array_values = array_values,
            array_value_ids_short = array_value_ids_short,
            dataframes = dataframes,
            dataframe_ids_short = dataframe_ids_short
        )

        ### for now, just check if the items in featurelist are really in the data frame and if not, skip them

        df_pl12 = df_pl12.fillna(value = 0) # impute missig values by 0

        df_pl12_scaled = pd.DataFrame(StandardScaler().fit_transform(df_pl12)) 

        #perform PCA
        pca = PCA(n_components=2)

        principalComponents = pca.fit_transform(df_pl12_scaled)

        col = range(1,(principalComponents.shape[1]+1))
        col = ['pca'+ str(y) for y in col]
        principalDf = pd.DataFrame(data = principalComponents, columns = col)

        principalDf["t_scaled"] = t_scaled
        principalDf["raw_file"] = hdf5_file_names

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
                        color_continuous_scale="bluered", title = "PCA on all data",
                        hover_name="raw_file", hover_data=["pca1", "pca2"],
                            labels={
                            "pca1": label_x,
                            "pca2": label_y,
                            "t_scaled": "timestamp"
                        })
            
        fig12.update_layout(width = int(args.width_pca), height = int(args.height_pca))
        fig12.update_traces(marker=dict(size=20))
            
        # Table and plot with feature loadings (weights of the variables in the PCA)
        loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index = df_pl12.columns)
        loadings.insert(0, "length", np.sqrt(loadings["PC1"]**2+ loadings["PC2"]**2))
        loadings.insert(0, "variable", loadings.index)
        loadings.sort_values("length", ascending=False, inplace=True)
        fig12_loadings = px.scatter(loadings, x = "PC1", y = "PC2", title = "PCA loadings (all data)", 
            hover_name="variable", hover_data=["PC1", "PC2"],)
        fig12_loadings.update_layout(width = int(args.width_pca), height = int(args.height_pca))
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
            width=args.width_pca,
            height=args.height_pca,
            title="Empty Plot"
        )
        fig12_loadings = fig12
        loadings = pd.DataFrame(columns = ["variable", "length", "PC1", "PC2"])
    if fig_show:
        fig12.show()
    with open(output_path + os.sep + "fig12a_PCA_all.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig12))
    fig12.write_html(file = output_path + os.sep + "fig12a_PCA_all.html", auto_open = False)
    if fig_show: 
        fig12_loadings.show()
    with open(output_path + os.sep + "fig12b_Loadings_all.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig11_loadings))
    fig12_loadings.write_html(file = output_path + os.sep + "fig12b_Loadings_all.html", auto_open = False)
   
    ### save loadings as tables
    if args.output_table_type == "csv":
        loadings.to_csv(output_path + os.sep + "fig12c_table_loadings_raw.csv", index = False)
    if args.output_table_type == "tsv":   
        loadings.to_csv(output_path + os.sep + "fig12c_table_loadings_raw.tsv", index = False, sep = "\t")
    if args.output_table_type == "xlsx":
        loadings.to_excel(output_path + os.sep + "fig12c_table_loadings_raw.xlsx", index = False)   



#################################################################################################
    ### Fig 13: Ion Maps (one for each raw file)

    if not os.path.exists(output_path + os.sep + "fig13_MS1_map"):
        os.makedirs(output_path + os.sep + "fig13_MS1_map")

    for file in hdf5_file_names:
        df_MS1_map = dataframes[file]["MS1_map"]

        ### reduce the number of points to roughly 1 million
        if (len(df_MS1_map) > 1000000):
                samples = int(len(df_MS1_map) / 1000000)
                df_MS1_map2 = df_MS1_map.loc[range(0, len(df_MS1_map), samples),:]
        else: 
                df_MS1_map2 = df_MS1_map

        df_MS1_map2["log_intensity"] = np.log10(df_MS1_map2["intensity"])

        if args.RT_unit == "min":
            df_MS1_map2["retention_time"] = df_MS1_map2["retention_time"]/60

        fig,ax = plt.subplots(figsize=(15,6)) 
        points = ax.scatter(df_MS1_map2["retention_time"], df_MS1_map2["mz"], c=df_MS1_map2["log_intensity"], s=1, cmap="Blues")
        fig.colorbar(points, label = "log10_intensity")
        fig.set_figheight(int(args.width_ionmaps))
        fig.set_figwidth(int(args.height_ionmaps))
        
        if args.RT_unit == "sec":
            ax.set_xlabel("retention time (sec)")
        elif args.RT_unit == "min":
            ax.set_xlabel("retention time (min)")
        ax.set_ylabel("m/z")
        ax.set_title(file)
        if fig_show:
            fig.show()
        fig.savefig(output_path + os.sep + "fig13_MS1_map" + os.sep + "fig13_MS1_map_" + file + ".png")

        
        

        

################################################################################################
    ### Figure 14: Pump Pressure
    
  
    ### start with empty plot, that is overwritten if pump pressure data is available
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
    
    pump_df = []
    for file in hdf5_file_names:
        if ("Pump_Pressure" not in dataframes[file].keys()):
            # Skip, there is no data available for this hdf5 file
            continue
        
        df_tmp_long = dataframes[file]["Pump_Pressure"]
        # use not more than roughly 10000 data points. If data has more than 10000 data points, take every nth data point        
        if df_tmp_long.shape[0] > 10000: 
            samples = int(df_tmp_long.shape[0] / 10000)
            df_tmp_long = df_tmp_long.iloc[range(0, df_tmp_long.shape[0], samples)]
        
        df_tmp_long = df_tmp_long.assign(filename=file)
        pump_df.append(df_tmp_long)    
           
 
        
    if (not pump_df == []):

        df_fig14_long = pd.concat(pump_df)
        ### x Axis data for pump pressure are in minutes, convert to seconds if necessary
        ### (this only holds for Thermo, for Bruker something is strange -> TODO)
        if args.RT_unit == "sec":
            df_fig14_long["pump_pressure_x_axis"] = df_fig14_long["pump_pressure_x_axis"]*60
        
        fig14 = px.line(df_fig14_long, x="pump_pressure_x_axis", y="pump_pressure_y_axis", color = "filename", title = "Pump Pressure")
        fig14.update_traces(line=dict(width=0.5))
        fig14.update_yaxes(exponentformat="E") 
        fig14.update_layout(height = int(args.height_barplots))
        if args.width_barplots > 0:
            fig14.update_layout(width = int(args.width_barplots))
        fig14.update_layout(yaxis_title = "Pump pressure")
        
        if args.RT_unit == "sec":
            fig14.update_layout(xaxis_title = "Time (sec)")
        elif args.RT_unit == "min":
            fig14.update_layout(xaxis_title = "Time (min)")
        
        
        fig14.write_html(file = output_path + os.sep + "fig14_Pump_pressure.html", auto_open = False)

    if fig_show:
        fig14.show()
    with open(output_path + os.sep + "fig14_Pump_pressure.plotly.json", "w") as json_file:
        json_file.write(plotly.io.to_json(fig14))
    fig14.write_html(file = output_path + os.sep + "fig14_Pump_pressure.html", auto_open = False)




################################################################################################
## Fig 15: all other headers by Thermo or Bruker

    if not os.path.exists(output_path + os.sep + "fig15_additional_headers"):
        os.makedirs(output_path + os.sep + "fig15_additional_headers")

    add_headers = []
    for file in hdf5_file_names:
        if ("Extracted_Headers" not in dataframes[file].keys()):
            continue
        else: 
            add_headers.extend(dataframes[file]["Extracted_Headers"].columns)
            
    if add_headers != []:
        ### remove duplicates
        add_headers = set(add_headers)
        add_headers = list(add_headers)# .sort()  
        add_headers.sort() ## sort alphabetically
        
        ### headers that define the time (x-axis)
        time_header = ""
        if "Time" in add_headers:
            time_header = "Time" # Bruker: Time 
        if "Scan_StartTime" in add_headers:
            time_header = "Scan_StartTime" # Thermo: Scan_StartTime

        if time_header == "":
            print("No time header found in the extracted headers, cannot plot additional headers!")
        else: 
    
            for header in add_headers:
                
                if header == time_header: # Time, will be needed as x-axis in all plots
                    continue
                if header == "MsMsType":  # MsMsType codes for DIA/DDA for example. Doesn't need to be plotted.
                    continue
                if header == "Scan_msLevel":  # MsMsType codes for MS1 or MS2 level. Doesn't need to be plotted.
                    continue
                
                ### extract data from compressed columns and put them into long format
                x = [] # x-axis time_header
                y = [] # y-axis additional header
                fn = [] # filename
                
                for file in hdf5_file_names:

                    if header not in dataframes[file]["Extracted_Headers"].columns:
                        # Skip, there is no data available for this hdf5 file
                        continue
                    
                    y_tmp = dataframes[file]["Extracted_Headers"][header].values
                    x_tmp = dataframes[file]["Extracted_Headers"][time_header].values
                    
                    display_header = header
                    # Ion injection time and lock mass correction should be filtered to only contain values for MS1 spectra
                    if header in ["EXTRA_Ion Injection Time (ms)", "EXTRA_LM mz-Correction (ppm),LM Correction"]:
                        msmsLevel = dataframes[file]["Extracted_Headers"]["Scan_msLevel"].values
                        y_tmp = y_tmp[msmsLevel == 1]
                        x_tmp = x_tmp[msmsLevel == 1]    
                        display_header = header + " (MS1 filtered)"
                    
                    x += [float(_x) for _x in x_tmp]
                    y += [float(_y) for _y in y_tmp]
                    fn += [file] * len(x_tmp)

                df_tmp = pd.DataFrame({
                    "filename": fn,
                    "x": x,
                    "y": y
                })
                
                ### transform time to minutes if necessary
                if args.RT_unit == "min":
                    df_tmp["x"] = df_tmp["x"]/60
                
            
                if not df_tmp.empty:
                    #df_tmp = df_tmp.sort_values(by = ["filename", "x"], ascending=True)  
                    fig15 = px.line(df_tmp, x="x", y="y", color = "filename", title = display_header)
                    fig15.update_traces(line=dict(width=0.5))
                    fig15.update_yaxes(exponentformat="E") 
                    fig15.update_layout(height = int(args.height_barplots))
                    if args.width_barplots > 0:
                        fig15.update_layout(width = int(args.width_barplots))
                    fig15.update_layout(yaxis_title = display_header)
                    if args.RT_unit == "sec":
                        fig15.update_layout(xaxis_title = "Time (sec)")
                    elif args.RT_unit == "min":
                        fig15.update_layout(xaxis_title = "Time (min)")
                    
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

                if fig_show:
                    fig15.show()
                with open(output_path + os.sep + "fig15_additional_headers" + os.sep + "{}.plotly.json".format(re.sub('\W+','', display_header)), "w") as json_file:
                    json_file.write(plotly.io.to_json(fig15))
                fig15.write_html(file = output_path + os.sep + "fig15_additional_headers" + os.sep + "{}.html".format(re.sub('\W+','', display_header)), auto_open = False)

################################################################################################
## Fig 16: all other headers by Thermo or Bruker      
        
    if not os.path.exists(output_path + os.sep + "fig16_BRUKER_calibrants"):
        os.makedirs(output_path + os.sep + "fig16_BRUKER_calibrants")
        
    df_calibrants = pd.DataFrame()
    for file in hdf5_file_names:
        if ("Calibrants" not in dataframes[file].keys()):
            continue
        else: 
            df_calibrants_tmp = dataframes[file]["Calibrants"]
            df_calibrants_tmp["filename"] = file
            df_calibrants = pd.concat([df_calibrants, df_calibrants_tmp])
       
    if not df_calibrants.empty:     
        calibrants = df_calibrants[["calibrant_mz", "calibrant_mobility"]]
        calibrants = calibrants.drop_duplicates()
        
        ## loop over calibrants and plot them
        i = 0
        for index, row in calibrants.iterrows():    
            i = i + 1
            
            df_tmp = df_calibrants[df_calibrants["calibrant_mz"] == row["calibrant_mz"]]
            df_tmp = df_tmp[df_tmp["calibrant_mobility"] == row["calibrant_mobility"]]
            print(df_tmp)
            
            mz_tmp = row["calibrant_mz"]
            mobility_tmp = row["calibrant_mobility"]
            
            ### rt is given in milliseconds here, convert to seconds or minutes
            if args.RT_unit == "sec":
                df_tmp["observed_calibrant_rt"] = df_tmp["observed_calibrant_rt"]/1000
            elif args.RT_unit == "min":
                df_tmp["observed_calibrant_rt"] = df_tmp["observed_calibrant_rt"]/60000


            title_tmp = "Calibrant " + str(i) + " m/z: " + str(mz_tmp) + ", ion mobility: " + str(mobility_tmp)

            fig16a = px.line(df_tmp, x="observed_calibrant_rt", y="observed_calibrant_mz", color = "filename", title = title_tmp)
            fig16a.update_traces(line=dict(width=0.5))
            fig16a.add_hline(y=mz_tmp)
            fig16a.update_layout(height = int(args.height_barplots))
            if args.width_barplots > 0:
                fig16a.update_layout(width = int(args.width_barplots))
            if args.RT_unit == "sec":
                fig16a.update_layout(xaxis_title = "Time (sec)")
            elif args.RT_unit == "min":
                fig16a.update_layout(xaxis_title = "Time (min)")
            with open(output_path + os.sep + "fig16_BRUKER_calibrants" + os.sep + "fig16a_Calibrant_mz_" + str(i) + ".plotly.json", "w") as json_file:
                json_file.write(plotly.io.to_json(fig16a))
            fig16a.write_html(file = output_path + os.sep + "fig16_BRUKER_calibrants" + os.sep + "fig16a_Calibrant_mz_" + str(i) + ".html", auto_open = False)

            
            fig16b = px.line(df_tmp, x="observed_calibrant_rt", y="observed_calibrant_mobility", color = "filename", title = title_tmp)
            fig16b.update_traces(line=dict(width=0.5))
            fig16b.add_hline(y=mobility_tmp)
            fig16b.update_layout(height = int(args.height_barplots))
            if args.width_barplots > 0:
                fig16b.update_layout(width = int(args.width_barplots))
            if args.RT_unit == "sec":
                fig16b.update_layout(xaxis_title = "Time (sec)")
            elif args.RT_unit == "min":
                fig16b.update_layout(xaxis_title = "Time (min)")
            with open(output_path + os.sep + "fig16_BRUKER_calibrants" + os.sep + "fig16b_Calibrant_ionmobility" + str(i) + ".plotly.json", "w") as json_file:
                json_file.write(plotly.io.to_json(fig16b))
            fig16b.write_html(file = output_path + os.sep + "fig16_BRUKER_calibrants" + os.sep + "fig16b_Calibrant_ionmobility" + str(i) + ".html", auto_open = False)

                
                

            
            
    


   
# %%
