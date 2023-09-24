#!/bin/env python

import sys
import os
import sqlite3
import argparse
from alphatims.bruker import TimsTOF
import alphatims.plotting
import zlib
import pickle
import base64
import csv

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d_folder", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_csv", help="The Output statistics CSV")
    parser.add_argument("-headers_to_parse", "-htp", help="The Headers to parse. Can be applied multiple times", action="append", default=[
        "Vacuum_CurrentFore",
        "Vacuum_Extra4thGauge",
        "Vacuum_CurrentHigh",
        "Vacuum_CurrentFunnel",
        "Digitizer_CurrentTemp",
        "TOF_DeviceTempCurrentValue1",
        "TOF_DeviceTempCurrentValue2"
    ])

    parser.add_argument("-headers_to_parse_column_name", "-cn", help="The Headers column name which was parsed (in same order). Can be applied multiple times", action="append", default=[
        "BRUKER_Vacuum_CurrentFore_pickle_zlib",
        "BRUKER_Vacuum_Extra4thGauge_pickle_zlib",
        "BRUKER_Vacuum_CurrentHigh_pickle_zlib",
        "BRUKER_Vacuum_CurrentFunnel_pickle_zlib",
        "BRUKER_Digitizer_CurrentTemp_pickle_zlib",
        "BRUKER_TOF_DeviceTempCurrentValue1_pickle_zlib",
        "BRUKER_TOF_DeviceTempCurrentValue2_pickle_zlib"
    ])
    return parser.parse_args()

if __name__ == "__main__":

    args = argparse_setup()

    args.d_folder = "/home/luxii/Documents/RAWs_BRUKER/TIM0000136_GA6_1_353.d"
    

    con = sqlite3.connect(args.d_folder + os.sep + "analysis.tdf")
    cur = con.cursor()

    # Get all property definitions
    res = cur.execute(
        "SELECT ID, PermanentName from PropertyDefinitions"
    )
    properties = res.fetchall()
    property_names = [x[1] for x in properties]

    # Filter to only needed ones
    p_index = []
    p_name = []
    p_col_name = []
    for n, col_n in zip(args.headers_to_parse, args.headers_to_parse_column_name):
        try:
            p_index.append(properties[property_names.index(n)][0])
            p_name.append(n)
            p_col_name.append(col_n)
        except:
            print("WARNING: Property '{}' not found!".format(n))


    # Extract data for each frame:
    data_dict = dict()
    for idx, name, col_name in zip(p_index, p_name, p_col_name):
        res = cur.execute(
            "SELECT Frame, Value from Properties WHERE  Property = {}".format(idx)
        )
        metadata = res.fetchall()
        data_dict[col_name] = [x[1] for x in sorted(metadata, key=lambda x: x[0])]

    # Special CASE: Get MS/MS-Type, Pressure and Measure/RetentionTime
    res = cur.execute(
        "SELECT Id, Time, MsMsType, Pressure from Frames"
    )
    frame_data = res.fetchall()
    sorted_frame_data = sorted(frame_data, key=lambda x: x[0])
    data_dict["BRUKER_Time_pickle_zlib"] = [x[1] for x in sorted_frame_data]
    data_dict["BRUKER_MsMsType_pickle_zlib"] = [x[2] for x in sorted_frame_data]
    data_dict["BRUKER_Pressure_pickle_zlib"] = [x[3] for x in sorted_frame_data]

    res.close()
    cur.close()
    con.close()


    # SPECIAL Case about PumpPreasure and so on
    # NC_Pump_Pressure
    con = sqlite3.connect(args.d_folder + os.sep + "chromatography-data.sqlite")
    cur = con.cursor()
    res = cur.execute(
        "SELECT Description, Id from TraceSources WHERE Description = 'NC_Pump_Pressure'"
    )
    traces = res.fetchall()

    if len(traces) != 0:
        trace_id = traces[0][1]
        res = cur.execute(
            "SELECT Times, Intensities from TraceChunks WHERE Trace = {}".format(trace_id)
        )
        trace_chunks = res.fetchall() 

        # Extract from binary data
        times = []
        data = []
        for t_chunk, d_chunk in trace_chunks:
            for t_off in range(0, len(t_chunk), 8):
                times.append(
                    int.from_bytes(t_chunk[t_off:t_off+8], byteorder="little")
                )
            for t_off in range(0, len(d_chunk), 4):
                data.append(
                    int.from_bytes(d_chunk[t_off:t_off+4], byteorder="little")
                )

        # TODO it is not clear which format the timestamps and the data has. All we know is that the data is provided in little endian
        data_dict["BRUKER_pump_preasure_bar_y_axis"] = times
        data_dict["BRUKER_pump_preasure_bar_x_axis"] = data

    # zlib and pickle all the data
    for k, v in data_dict.items():
        data_dict[k] = base64.b64encode(zlib.compress(pickle.dumps(v), level=9)).decode("utf-8")



    # Write Output
    with open(args.out_csv, "w") as csv_out:

        csv_writer = csv.DictWriter(csv_out, fieldnames=data_dict.keys())
        csv_writer.writeheader()
        csv_writer.writerow(data_dict)
