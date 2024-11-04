#!/usr/bin/env python

import argparse
import os
import sqlite3
from typing import Any, List

import h5py
import numpy as np


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d_folder", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_hdf5", help="The Output statistics HDF5")
    parser.add_argument("-headers_to_parse", "-htp", help="The Headers to parse. Can be applied multiple times", action="append", default=[
        "Vacuum_CurrentFore",
        "Vacuum_Extra4thGauge",
        "Vacuum_CurrentHigh",
        "Vacuum_CurrentFunnel",
        "Digitizer_CurrentTemp",
        "TOF_DeviceTempCurrentValue1",
        "TOF_DeviceTempCurrentValue2"
    ])

    return parser.parse_args()


def add_entry_to_hdf5(
    f, qc_acc: str, qc_short_name: str, qc_name: str, qc_description: str, 
    value, value_shape: tuple, value_type: str, 
    unit_accession: str=None, unit_name: str=None,
    ): 
    """ Adds an entry into the hdf5 file """
    key = "|".join([qc_acc, qc_short_name])  # ACCESSION|SHORT_DESC
    if value_type in ("str", h5py.string_dtype()):
        ds = f.create_dataset(key, shape=value_shape, dtype=h5py.string_dtype(), compression="gzip")
        ds[:] = value
    else:
        f.create_dataset(key, value_shape, dtype=value_type, compression="gzip")
        f[key].write_direct(np.array(value, dtype=value_type))
    
    f[key].attrs["qc_short_name"] = qc_short_name
    f[key].attrs["qc_name"] = qc_name
    f[key].attrs["qc_description"] = qc_description
    f[key].attrs["unit_accession"] = unit_accession
    f[key].attrs["unit_name"] = unit_name


def add_table_in_hdf5(
    f, qc_acc: str, qc_short_name: str, qc_name: str, qc_description: str, 
    column_names: List[str], column_data: List[List[Any]], column_types: List[str]
    ): 
    """Adds a table in groups"""

    key = "|".join([qc_acc, qc_short_name])  # ACCESSION|SHORT_DESC
    table_group = f.create_group(key)
    table_group.attrs["qc_short_name"] = qc_short_name
    table_group.attrs["qc_name"] = qc_name
    table_group.attrs["qc_description"] = qc_description
    table_group.attrs["column_order"] = "|".join(column_names)
    
    for n, d, t in zip (column_names, column_data, column_types): 
        if t in ("str", h5py.string_dtype()):
            ds = table_group.create_dataset(n, shape=len(d), dtype=h5py.string_dtype(), compression="gzip")
            ds[:] = d
        else:
            table_group.create_dataset(n, (len(d),), dtype=t, compression="gzip")
            table_group[n].write_direct(np.array(d, dtype=t))


if __name__ == "__main__":

    args = argparse_setup()
    con = sqlite3.connect(args.d_folder + os.sep + "analysis.tdf")
    cur = con.cursor()

    # Get all property definitions
    res = cur.execute(
        "SELECT ID, PermanentName from PropertyDefinitions"
    )
    properties = res.fetchall()
    property_names = [x[1] for x in properties]

    print("Available Properties: ")
    for p in property_names:
        print("\t" + p)

    # Filter to only needed ones
    p_index = []
    p_name = []
    p_col_name = []
    for n in args.headers_to_parse:
        try:
            p_index.append(properties[property_names.index(n)][0])
            p_name.append(n)
        except:
            print("WARNING: Property '{}' not found!".format(n))


    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:

        # Extract data for each frame:
        data_dict = dict()
        for idx, name in zip(p_index, p_name):
            res = cur.execute(
                "SELECT Frame, Value from Properties WHERE  Property = {}".format(idx)
            )
            metadata = res.fetchall()

            data_dict[name] = [x[1]if x[1] is not None else np.nan for x in sorted(metadata, key=lambda x: x[0])]

        # Special CASE: Get MS/MS-Type, Pressure and Measure/RetentionTime
        res = cur.execute(
            "SELECT Id, Time, MsMsType, Pressure from Frames"
        )
        frame_data = res.fetchall()
        sorted_frame_data = sorted(frame_data, key=lambda x: x[0])

        column_name = list(data_dict.keys()) + ["Time", "MsMsType", "Pressure"]
        column_data = [data_dict[x] for x in column_name[:-3]] + [[x[1] for x in sorted_frame_data], [x[2] for x in sorted_frame_data], [x[3] for x in sorted_frame_data]]
        column_type = ["float64"]*len(column_name)
        add_table_in_hdf5(
            out_h5, "BRUKER", "Extracted_Headers", "The extracted Bruker headers, which have been specified prior.", 
            "This table can contain various columns, like 'Temperature' and more. Depending "
            "on the input Bruker-file a column may be present in this table.",
            column_name, column_data, column_type
        )

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
        else:
            times, data = [np.nan], [np.nan]

        # TODO it is not clear which format the timestamps and the data has. All we know is that the data is provided in little endian
        column_name = ["pump_pressure_x_axis", "pump_pressure_y_axis"]
        column_data = [times, data]
        column_type = ["float64", "float64"]
        add_table_in_hdf5(
            out_h5, "LOCAL:21", "Pump_Pressure", "Pump Pressure",
            "The pump pressure as a table, containing the coordinates.",
            column_name, column_data, column_type
        )
