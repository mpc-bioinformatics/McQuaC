#!/bin/env python

import sys
import argparse
import pyopenms
import numpy as np
import math
from collections import defaultdict
import datetime
import time
import csv
import zlib
import pickle
import base64

from fisher_py import RawFile

from fisher_py.raw_file_reader import RawFileReaderAdapter, RawFileAccess
from fisher_py.data.business import GenericDataTypes, ChromatogramTraceSettings, TraceType, ChromatogramSignal, SpectrumPacketType, Scan
from fisher_py.data.filter_enums import MsOrderType
from fisher_py.data import Device, ToleranceUnits
from fisher_py.mass_precision_estimator import PrecisionEstimate


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-raw", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_csv", help="The Output statistics CSV")
    parser.add_argument("-headers_to_parse", "-htp", help="The Headers to parse. Can be applied multiple times", action="append", default=[
        "Ion Injection Time",
        "Number of Lock Masses",
        "Lock Mass #1",
        "Lock Mass #2",
        "Lock Mass #3",
        "LM Search Window",
        "LM Search Window",
        "Number of LM Found",
        "Last Locking",
        "LM m/z-Correction,LM Correction",
    ])
    parser.add_argument("-headers_to_parse_column_name", "-cn", help="The Headers column name which was parsed (in same order). Can be applied multiple times", action="append", default=[
        "THERMO_Ion_Injection_Time_pickle_zlib",
        "THERMO_Number_of_Lock_Masses_pickle_zlib",
        "THERMO_Lock_Mass_1_pickle_zlib",
        "THERMO_Lock_Mass_2_pickle_zlib",
        "THERMO_Lock_Mass_3_pickle_zlib",
        "THERMO_LM_Search_Window_pickle_zlib",
        "THERMO_LM_Search_Window_pickle_zlib",
        "THERMO_Number_of_LM_Found_pickle_zlib",
        "THERMO_Last_Locking_pickle_zlib",
        "THERMO_LM_m_z_Correction_pickle_zlib",
    ])
    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_setup()
    data_dict = defaultdict(lambda: list())

    raw_file = RawFileReaderAdapter.file_factory(args.raw)
    raw_file.select_instrument(Device.MS, 1)  # Selecting the MS
    all_headers =  raw_file.get_trailer_extra_header_information()

    statistics_to_retrieve = []
    print("RAW-file: {}".format(args.raw))
    for idx, h in enumerate(all_headers):

        print("Field {},\t\tLabel: {}".format(idx, h.label))
        for hp, cn in zip(args.headers_to_parse, args.headers_to_parse_column_name):
            for hp_part in hp.split(","):
                if h.label.startswith(hp_part):
                    # Retrievable in this RAW-file
                    statistics_to_retrieve.append((idx, hp, cn))
                    break


    # Retrieve all the information 
    first_scan_number = raw_file.run_header_ex.first_spectrum
    last_scan_number = raw_file.run_header_ex.last_spectrum
    for scan in range(first_scan_number, last_scan_number + 1):
        scan_values = raw_file.get_trailer_extra_information(scan).values
        scan_statistics = raw_file.get_scan_stats_for_scan_number(scan)
        start_time_of_scan = scan_statistics.start_time
        data_dict["THERMO_Scan_StartTime_zlib"].append(start_time_of_scan)
        scan_filter = raw_file.get_filter_for_scan_number(scan)
        data_dict["THERMO_Scan_msLevel_zlib"].append(scan_filter.ms_order.value)

        # Get Info of filtered statistics we want to track
        for idx, hp, cn in statistics_to_retrieve:
            data_dict[cn].append(scan_values[idx])


    # Get all the information from the instruments (in FreeStyle under Devices)
    num_devices = raw_file.get_instrument_count_of_type(Device.Analog)  # Get the number of devices
    settings = ChromatogramTraceSettings(TraceType.Analog1)
    
    print("\nAnalog Devices:")
    for i in range(1, num_devices+1):
        # Iterate over each device
        raw_file.select_instrument(Device.Analog, i)
        label = raw_file.get_instrument_data().axis_label_y
        
        print("Info: {} (Label: {})".format(label, raw_file.get_instrument_data().channel_labels[0]))

        # Check if it is the pump_preasure, which we want to extract
        if label in ("Pump_Pressure bar", "NC_Pump_Pressure bar"):
            # This works for QexHF, QeXI and FLI (NC_Pump_Pressure) and EX, EXI and EXII (Pump_Preasure)
            data = raw_file.get_chromatogram_data([settings], 1, -1)
            data_dict["THERMO_pump_preasure_bar_y_axis"] = list(data.intensities_array[0])
            data_dict["THERMO_pump_preasure_bar_x_axis"] = list(data.positions_array[0])
   
    if "pump_preasure_bar_x_axis" not in data_dict:
        # This is currently the case for PROETD and OEI. These are not captured, hence None
        data_dict["THERMO_pump_preasure_bar_x_axis"] = None
        data_dict["THERMO_pump_preasure_bar_y_axis"] = None

    # Close raw-file
    raw_file.dispose()  

    # zlib and pickle all the data
    for k, v in data_dict.items():
        data_dict[k] = base64.b64encode(zlib.compress(pickle.dumps(v), level=9)).decode("utf-8")

    # Write Output
    with open(args.out_csv, "w") as csv_out:

        csv_writer = csv.DictWriter(csv_out, fieldnames=data_dict.keys())
        csv_writer.writeheader()
        csv_writer.writerow(data_dict)
