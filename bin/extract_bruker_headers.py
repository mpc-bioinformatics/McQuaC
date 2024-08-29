#!/usr/bin/env python

import os
import sqlite3
import argparse
import numpy as np
import h5py


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


def add_entry_to_hdf5(f, key: str, value, array_shape: tuple, data_type: str, unit: str, description: str, compression: str=None):
    """ Adds an entry into the hdf5 file """
    out_h5.create_dataset(key, array_shape, dtype=data_type, compression=compression)
    out_h5[key].attrs["unit"] = unit
    out_h5[key].attrs["Description"] = description
    out_h5[key].write_direct(np.array(value, dtype=data_type))


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
            p_col_name.append("BRUKER_" + n)
        except:
            print("WARNING: Property '{}' not found!".format(n))


    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'a') as out_h5:

        # Extract data for each frame:
        for idx, name, col_name in zip(p_index, p_name, p_col_name):
            res = cur.execute(
                "SELECT Frame, Value from Properties WHERE  Property = {}".format(idx)
            )
            metadata = res.fetchall()

            entries = [x[1]if x[1] is not None else np.nan for x in sorted(metadata, key=lambda x: x[0])]
            add_entry_to_hdf5(
                out_h5, col_name, entries, (len(entries),), "float64", "refer to description", 
                description=name
            )

        # Special CASE: Get MS/MS-Type, Pressure and Measure/RetentionTime
        res = cur.execute(
            "SELECT Id, Time, MsMsType, Pressure from Frames"
        )
        frame_data = res.fetchall()
        sorted_frame_data = sorted(frame_data, key=lambda x: x[0])

        len_entries = len([x[1] for x in sorted_frame_data])
        add_entry_to_hdf5(
            out_h5, "BRUKER_Time", [x[1] for x in sorted_frame_data], (len_entries,), "float64", "seconds", 
            description="Timepoints for BRUKER RAW-arrays in seconds"
        )
        add_entry_to_hdf5(
            out_h5, "BRUKER_MsMsType", [x[2] for x in sorted_frame_data], (len_entries,), "float64", "none", 
            description="MsMsType ( 8 --> ddaPASEF, 9 --> diaPASSEF)"
        )
        add_entry_to_hdf5(
            out_h5, "BRUKER_Pressure", [x[3] for x in sorted_frame_data], (len_entries,), "float64", "refer to description", 
            description="Pressure"
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
        add_entry_to_hdf5(
            out_h5, "BRUKER_pump_pressure_bar_x_axis", times, (len(times),), "float64", "unknown", 
            description="Pump Pressure x coordinate, NaN if not present", compression="gzip"
        )
        add_entry_to_hdf5(
            out_h5, "BRUKER_pump_pressure_bar_y_axis", data, (len(data),), "float64", "unknown", 
            description="Pump Pressure y coordinate, NaN if not present", compression="gzip"
        )
