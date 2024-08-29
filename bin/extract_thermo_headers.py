#!/usr/bin/env python

import argparse
import pyopenms
import numpy as np
import math
from collections import defaultdict
import datetime
import time
import h5py

from fisher_py import RawFile

from fisher_py.raw_file_reader import RawFileReaderAdapter, RawFileAccess
from fisher_py.data.business import GenericDataTypes, ChromatogramTraceSettings, TraceType, ChromatogramSignal, SpectrumPacketType, Scan
from fisher_py.data.filter_enums import MsOrderType
from fisher_py.data import Device, ToleranceUnits
from fisher_py.mass_precision_estimator import PrecisionEstimate


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


def add_entry_to_hdf5(f, key: str, value, array_shape: tuple, data_type: str, unit: str, description: str, compression: str=None):
    """ Adds an entry into the hdf5 file """
    out_h5.create_dataset(key, array_shape, dtype=data_type, compression=compression)
    out_h5[key].attrs["unit"] = unit
    out_h5[key].attrs["Description"] = description
    out_h5[key].write_direct(np.array(value, dtype=data_type))


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
    with h5py.File(args.out_hdf5, 'a') as out_h5:

        # Retrieve all the information 
        first_scan_number = raw_file.run_header_ex.first_spectrum
        last_scan_number = raw_file.run_header_ex.last_spectrum
        for scan in range(first_scan_number, last_scan_number + 1):
            scan_statistics = raw_file.get_scan_stats_for_scan_number(scan)
            start_time_of_scan = scan_statistics.start_time
            data_dict["THERMO_Scan_StartTime"].append(start_time_of_scan)
            scan_filter = raw_file.get_filter_for_scan_number(scan)
            data_dict["THERMO_Scan_msLevel"].append(scan_filter.ms_order.value)

            # Get Info of filtered statistics we want to track (log and extra data)
            log_scan_values = raw_file.get_status_log_for_retention_time(start_time_of_scan).values
            for idx, hp in log_statistics_to_retrieve:
                data_dict["THERMO_LOG_" + hp].append(log_scan_values[idx])
            extra_scan_values = raw_file.get_trailer_extra_information(scan).values
            for idx, hp in extra_statistics_to_retrieve:
                data_dict["THERMO_EXTRA_" + hp].append(extra_scan_values[idx])


        for key, val in data_dict.items():
            if "THERMO_LOG_" in key:
                desc = key[len("THERMO_LOG_"):] 
            elif "THERMO_EXTRA_" in key:
                desc = key[len("THERMO_EXTRA_"):]
            else:
                desc = key[len("THERMO_"):]
            add_entry_to_hdf5(
                out_h5, key.replace("/", ""), val, (len(val),), "float64", "refer to description", 
                description=desc, compression="gzip"
            )


        # Tune data is one dimensional, therefore single values
        tune_scan_values = raw_file.get_tune_data(0).values
        for idx, hp in tune_statistics_to_retrieve:
            add_entry_to_hdf5(
                out_h5, "THERMO_TUNE_" + hp, tune_scan_values[idx], (1,), "float64", "refer to description", 
                description=hp
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
                pump_pressure_dict["THERMO_pump_pressure_bar_y_axis"] = list(data.intensities_array[0])
                pump_pressure_dict["THERMO_pump_pressure_bar_x_axis"] = list(data.positions_array[0])
    
        if "THERMO_pump_pressure_bar_x_axis" not in pump_pressure_dict:
            # This is currently the case for PROETD and OEI. These are not captured, hence None
            pump_pressure_dict["THERMO_pump_pressure_bar_x_axis"] = [np.nan]
            pump_pressure_dict["THERMO_pump_pressure_bar_y_axis"] = [np.nan]

        entries_in_pp = len(pump_pressure_dict["THERMO_pump_pressure_bar_x_axis"])
        add_entry_to_hdf5(
            out_h5, "THERMO_pump_pressure_bar_x_axis", pump_pressure_dict["THERMO_pump_pressure_bar_x_axis"], (entries_in_pp,), "float64", "minutes", 
            description="Times when pump pressure was measured (is NaN if nothing was recorded)", compression="gzip"
        )
        add_entry_to_hdf5(
            out_h5, "THERMO_pump_pressure_bar_y_axis", pump_pressure_dict["THERMO_pump_pressure_bar_x_axis"], (entries_in_pp,), "float64", "refer to description", 
            description="The pump pressure over time  (is NaN if nothing was recorded)", compression="gzip"
        )

    # Close raw-file
    raw_file.dispose()  
