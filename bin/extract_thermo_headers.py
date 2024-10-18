#!/usr/bin/env python

import argparse
import datetime
import math
import time
from collections import defaultdict
from typing import Any, List

import h5py
import numpy as np
import pyopenms
from fisher_py import RawFile
from fisher_py.data import Device, ToleranceUnits
from fisher_py.data.business import (ChromatogramSignal,
                                     ChromatogramTraceSettings,
                                     GenericDataTypes, Scan,
                                     SpectrumPacketType, TraceType)
from fisher_py.data.filter_enums import MsOrderType
from fisher_py.mass_precision_estimator import PrecisionEstimate
from fisher_py.raw_file_reader import RawFileAccess, RawFileReaderAdapter


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-raw", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_hdf5", help="The Output statistics HDF5")
    parser.add_argument("-extra_headers_to_parse", "-ehtp", help="The Headers to parse. Can be applied multiple times", action="append", 
        default=[
            "Ion Injection Time (ms)",
            "Number of Lock Masses",
            "Lock Mass #1 (m/z)",
            "Lock Mass #2 (m/z)",
            "Lock Mass #3 (m/z)",
            "LM Search Window (ppm)",
            "LM Search Window (mmu)",
            "Last Locking (sec)",
            "LM m/z-Correction (ppm),LM Correction",
    ])
    parser.add_argument("-tune_headers_to_parse", "-thtp", help="The Headers to parse. Can be applied multiple times. NOTE: These entries are one dimensional.", action="append", 
        default=[
            "Ion Transfer Tube Temperature (+ or +-)",
            "Ion Transfer Tube Temperature (-)",
            "Vaporizer Temp. (+ or +-)",
            "Vaporizer Temp. (-)",
    ])
    parser.add_argument("-log_headers_to_parse", "-lhtp", help="The Headers to parse. Can be applied multiple times", action="append", 
        default=[
            "Ion Transfer Tube Temperature (°C)",
            "Vaporizer Temperature (°C)",
            "Ambient temp. (°C)",
            "Orbitrap block temp. (°C)",
            "Detector temp. (°C)",
            "Ion tr. tube temp. (°C)",
            "Vaporizer temp. (°C)",
            "CPU core temp. (°C)",
            "PCB temp. top (°C)",
            "PCB temp. center (°C)",
            "PCB temp. bottom (°C)",
            "Inner Electrode Temp. (°C)",
            "Outer Electrode 1 Temp. (°C)",
            "Outer Electrode 2 Temp. (°C)",
            "MCB Ambient Temp. (°C)",
            "MCB PCB Temp. (°C)",
            "TX PCB temperature (°C)",
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


def get_headers_to_parse(headers, headers_from_raw):
    """ Helper function to map which header from which index should be retrieved"""
    statistics_to_retrieve = []
    for idx, h in enumerate(headers_from_raw):
        print("Field {},\t\tLabel: {}".format(idx, h.label))
        for hp in headers:
            for hp_part in hp.split(","):
                if h.label.startswith(hp_part):
                    # Retrievable in this RAW-file
                    statistics_to_retrieve.append((idx, hp))
                    break

    return statistics_to_retrieve


if __name__ == "__main__":
    args = argparse_setup()

    data_dict = defaultdict(lambda: list())
    raw_file = RawFileReaderAdapter.file_factory(args.raw)
    raw_file.select_instrument(Device.MS, 1)  # Selecting the MS

    print("Retrieveing log headers")
    log_statistics_to_retrieve = get_headers_to_parse(
        args.log_headers_to_parse,
        raw_file.get_status_log_header_information()
    )
    print("Retrieveing tune headers")
    tune_statistics_to_retrieve = get_headers_to_parse(
        args.tune_headers_to_parse,
        raw_file.get_tune_data_header_information()
    )
    print("Retrieveing extra headers")
    extra_statistics_to_retrieve = get_headers_to_parse(
        args.extra_headers_to_parse,
        raw_file.get_trailer_extra_header_information()
    )

    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:

        # Retrieve all the information 
        first_scan_number = raw_file.run_header_ex.first_spectrum
        last_scan_number = raw_file.run_header_ex.last_spectrum
        for scan in range(first_scan_number, last_scan_number + 1):
            scan_statistics = raw_file.get_scan_stats_for_scan_number(scan)
            start_time_of_scan = scan_statistics.start_time
            data_dict["Scan_StartTime"].append(start_time_of_scan)
            scan_filter = raw_file.get_filter_for_scan_number(scan)
            data_dict["Scan_msLevel"].append(scan_filter.ms_order.value)

            # Get Info of filtered statistics we want to track (log and extra data)
            log_scan_values = raw_file.get_status_log_for_retention_time(start_time_of_scan).values
            for idx, hp in log_statistics_to_retrieve:
                data_dict["LOG_" + hp].append(log_scan_values[idx])
            extra_scan_values = raw_file.get_trailer_extra_information(scan).values
            for idx, hp in extra_statistics_to_retrieve:
                data_dict["EXTRA_" + hp].append(extra_scan_values[idx])


        if len(list(data_dict.keys())) != 0:
            column_name = list(data_dict.keys())
            column_data = [data_dict[x] for x in column_name]
            column_type = ["float64"]*len(column_name)
            add_table_in_hdf5(
                out_h5, "THERMO", "Extracted_Headers", "The extracted Thermo headers, which have been specified.", 
                "This table can contain various columns, ranging from 'Temperature', 'Lock Masses' and more. Depending "
                "on the input RAW-file a column may be present in this table.",
                [x.replace("/", "") for x in column_name], column_data, column_type
            )

        # Tune data is one dimensional, therefore single values
        tune_scan_values = raw_file.get_tune_data(0).values
        if len(tune_statistics_to_retrieve) != 0:
            tune_dict = dict()
            for idx, hp in tune_statistics_to_retrieve:
                tune_dict["TUNE_" + hp], tune_scan_values[idx]
            column_name = list(tune_dict.keys())
            column_data = [tune_dict[x] for x in column_name]
            column_type = ["float64"]*len(column_name)
            add_table_in_hdf5(
                out_h5, "THERMO_LOG", "Extracted_Log_Headers", "The extracted Thermo log headers, which have been specified.", 
                "This table contains various columns of length one, like 'Vaporizer Temperature'. Depending "
                "on the input RAW-file a column may be present in this table.",
                [x.replace("/", "") for x in column_name], column_data, column_type
            )

        # Get all the information from the instruments (in FreeStyle under Devices)
        num_devices = raw_file.get_instrument_count_of_type(Device.Analog)  # Get the number of devices
        settings = ChromatogramTraceSettings(TraceType.Analog1)
        pump_pressure_dict = defaultdict(lambda: list())

        
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
                pump_pressure_dict["pump_pressure_bar_y_axis"] = list(data.intensities_array[0])
                pump_pressure_dict["pump_pressure_bar_x_axis"] = list(data.positions_array[0])

        if "pump_pressure_bar_y_axis" not in pump_pressure_dict:
            # This is currently the case for PROETD and OEI. These are not captured, hence None
            pump_pressure_dict["pump_pressure_bar_x_axis"] = [np.nan]
            pump_pressure_dict["pump_pressure_bar_y_axis"] = [np.nan]

        column_name = ["pump_pressure_x_axis", "pump_pressure_y_axis"]
        column_data = [pump_pressure_dict["pump_pressure_bar_x_axis"], pump_pressure_dict["pump_pressure_bar_y_axis"]]
        column_type = ["float64", "float64"]
        add_table_in_hdf5(
            out_h5, "LOCAL:21", "Pump_Pressure", "Pump Pressure",
            "The pump pressure as a table, containing the coordinates.",
            column_name, column_data, column_type
        )

    # Close raw-file
    raw_file.dispose()  
