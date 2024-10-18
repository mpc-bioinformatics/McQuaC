#!/usr/bin/env python

import sys
import argparse
# import pyopenms
import numpy as np
import math
from collections import defaultdict
import datetime
import time
import h5py
from typing import List, Any


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-mzml", help="mzML file to extract data from")
    parser.add_argument("-out_hdf5", help="The Output HDF5 file")
    parser.add_argument("-base_peak_tic_up_to", default=105, type=float, help="Retrieve the Base peak Intensity Max and the Total Ion "
        "Current from minute 0 up to minute X. Defaults to 105 minutes")
    parser.add_argument("-filter_threshold", default=0.00001, type=float, help="Threshold for the MS1 peaks, to be"
        " included in the output file. Defaults to 0.00001 (0.001%) of the highest overall MS1 peak. Values lower"
        " will be disregarded.")
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
    




# Total Ion Current Calculation taken from:
# https://pyopenms.readthedocs.io/en/latest/first_steps.html?highlight=TIC#total-ion-current-calculation
# This Could be the columns "accumulated-MS1_TIC" and "accumulated-MS2_TIC"
def get_accumulated_TIC(exp, mslevel):
    return sum(
        [
            sum(s.get_peaks()[1]) # Simply sum all the peaks up
            for s in exp if s.getMSLevel() == mslevel # Filter for only summing up of "a" MSLevel
        ]
    )


if __name__ == "__main__":
    # args = argparse_setup()

    # # Load MZML
    # exp = pyopenms.MSExperiment()
    # pyopenms.MzMLFile().load(args.mzml, exp)

    # Open HDF5 file in write mode
    with h5py.File("test.h5", 'w') as out_h5:

        # Get Date and Time of Measurement (timestamp)
        hh_mm_ss_str = exp.getDateTime().getTime()
        yy_mm_dd_str = exp.getDateTime().getDate() # could also be yy_mm_dd
        try:
            add_entry_to_hdf5(
                out_h5, "timestamp", time.mktime(datetime.datetime.strptime(yy_mm_dd_str + "|" + hh_mm_ss_str , "%Y-%m-%d|%H:%M:%S").timetuple()),
                (1,), "int64", "UNIX timestamp",
                "Date and time of the measurement as a UNIX timestamp. It is set to -1 if no data is available."
            ) 
        except:
            add_entry_to_hdf5(
                out_h5, "timestamp", -1,
                (1,), "int64", "UNIX timestamp",
                "Date and time of the measurement as a UNIX timestamp. It is set to -1 if no data is available."
            ) 

        # Get accumulated TICs for ms1 and ms2 (accumulated-MS1_TIC)
        add_entry_to_hdf5(
            out_h5, "accumulated-MS1_TIC", get_accumulated_TIC(exp, 1), (1,), "float64", "none",
            "The overall sum over all peaks across all MS1 spectra."
        ) 
        add_entry_to_hdf5(
            out_h5, "accumulated-MS2_TIC", get_accumulated_TIC(exp, 2), (1,), "float64", "none",
            "The overall sum over all peaks across all MS2 spectra."
        ) 

        # Get Number of MS 1 and 2 Spectra and also the precursors of the MS2 by counted charge
        num_ms1_spectra = 0
        num_ms2_spectra = 0
        num_ms2_prec_charges = defaultdict(lambda: 0) 
        for spectrum in exp.getSpectra():
            if spectrum.getMSLevel() == 1:
                num_ms1_spectra += 1
            elif spectrum.getMSLevel() == 2:
                num_ms2_spectra +=1
                num_ms2_prec_charges[spectrum.getPrecursors()[0].getCharge()] += 1

        add_entry_to_hdf5(
            out_h5, "total_num_ms1", num_ms1_spectra, (1,), "int32", "none",
            "Number of measured MS1 spectra."
        ) 
        add_entry_to_hdf5(
            out_h5, "total_num_ms2", num_ms2_spectra, (1,), "int32", "none",
            "Number of measured MS2 spectra"
        ) 

        prec_unknown = 0 if 0 not in num_ms2_prec_charges else num_ms2_prec_charges[0] / num_ms2_spectra    
        add_entry_to_hdf5(
            out_h5, "MS2_PrecZ_Unknown", prec_unknown, (1,), "float64", "none",
            "Proportion of MS2 precursors with unknown charge state (may happen for older machines)."
        )
        for i in range(1, 6):
            add_entry_to_hdf5(
                out_h5, "MS2_PrecZ_" + str(i), num_ms2_prec_charges[i] / num_ms2_spectra, (1,), "float64", "none",
                "Proportion of MS2 precursors with charge " + str(1) + "."
            )
        
        precz_more = 0
        for key in num_ms2_prec_charges.keys():
            if key > 5:
                precz_more += num_ms2_prec_charges[key]
        add_entry_to_hdf5(
            out_h5, "MS2_PrecZ_more", precz_more / num_ms2_spectra, (1,), "float64", "none",
            "Proportion of MS2 precursors with charge 6 and more."
        )
        
        # Here we have the last spectrum in "spectrum". We get the Retention Time Duration from it (in seconds)
        rt_duration = spectrum.getRT() # In Seconds: How long the MS has acquired spectra 
        add_entry_to_hdf5(
            out_h5, "RT_duration", rt_duration, (1,), "float64", "seconds",
            "Scan start time of the last measured spectrum."
        )

        # Get Base_peak_Intensity total_ion_current and up to
        ms1_ms2_basepeaks = []
        ms1_ms2_rt = []
        ms1_ms2_tic = []
        for spectrum in exp.getSpectra():
            ms1_ms2_rt.append(spectrum.getRT())
            peaks = spectrum.get_peaks()[1]
            ms1_ms2_basepeaks.append(max(peaks) if len(peaks) != 0 else 0)
            ms1_ms2_tic.append(sum(spectrum.get_peaks()[1]))

        base_peak_intensity_max = max(ms1_ms2_basepeaks)
        total_ion_current_max = max(ms1_ms2_tic)

        add_entry_to_hdf5(
            out_h5, "Base_Peak_Intensity_Max", base_peak_intensity_max, (1,), "float64", "none",
            "The maximum base peak (highest peak in spectrum) across all MS1 and MS2 spectra."
        )
        add_entry_to_hdf5(
            out_h5, "Total_Ion_Current_Max", total_ion_current_max, (1,), "float64", "none",
            "The maximum of all TICs across MS1 and MS2 spectra."
        )

        # and up to 105 minutes
        break_on = args.base_peak_tic_up_to*60
        for i in range(len(ms1_ms2_rt)):
            if break_on < ms1_ms2_rt[i]:
                break
        base_peak_intensity_max_up_to_xm = max(ms1_ms2_basepeaks[:i-1])
        total_ion_current_max_up_to_xm = max(ms1_ms2_tic[:i-1])

        add_entry_to_hdf5(
            out_h5, "Base_Peak_Intensity_Max_Up_To_" + str(args.base_peak_tic_up_to), base_peak_intensity_max, (1,), "float64", "none",
            "The maximum base peak (highest peak in spectrum) across all MS1 and MS2 spectra up to " + str(args.base_peak_tic_up_to) + " minutes."
        )
        add_entry_to_hdf5(
            out_h5, "Total_Ion_Current_Max_Up_To_" + str(args.base_peak_tic_up_to), total_ion_current_max, (1,), "float64", "none",
            "The maximum of all TICs across MS1 and MS2 spectra up to " + str(args.base_peak_tic_up_to) + " minutes."
        )
        ###

        # Get MS1 and MS2 TICs
        ms1_tic = [] # This is the MS1-TIC-BLOB 
        ms1_rt = []  # This is the MS1-RT-BLOB
        ms1_map_rt = []  # This is the MS1-RT-BLOB (for heatmap purpose, data is transformed to 3 dimensions)
        ms1_map_intens = []  # This is the MS1-RT-BLOB (for heatmap purpose, data is transformed to 3 dimensions)
        ms1_map_mz = []  # This is the MS1-RT-BLOB (for heatmap purpose, data is transformed to 3 dimensions)
        ms2_tic = [] # This is the MS2-TIC-BLOB
        ms2_rt = []  # This is the MS2-RT-BLOB
        ms2_mz = []  # This is the MS2-MZ-BLOB (needed or the Ion-Map)
        BASE_PEAK_NOISE_THRESHOLD = base_peak_intensity_max*args.filter_threshold  
        for spectrum in exp.getSpectra():
            if spectrum.getMSLevel() == 1:
                rt = spectrum.getRT()
                ms1_rt.append(rt)
                ms1_tic.append(sum(spectrum.get_peaks()[1]))

                # Filter by basepeak intensity 
                mz, intens = spectrum.get_peaks()
                for m, i in zip(mz, intens):
                    if i >= BASE_PEAK_NOISE_THRESHOLD:
                        ms1_map_mz.append(m)
                        ms1_map_intens.append(i)
                        ms1_map_rt.append(rt)
            elif spectrum.getMSLevel() == 2:
                ms2_rt.append(spectrum.getRT())
                ms2_tic.append(sum(spectrum.get_peaks()[1]))
                precs = spectrum.getPrecursors()
                if len(precs) != 1:
                    raise Exception("Unexpected number of precursors in a MS2-Spectrum")
                else:
                    ms2_mz.append(precs[0].getMZ())

        # Adding TIC and RT for TIC Plots
        add_entry_to_hdf5(
            out_h5, "ms1_tic_array", ms1_tic, (len(ms1_tic),), "float64", "none",
            "RAW-Array of the TICs of MS1 spectra (to plot tic overlay plot in combination with 'ms1_rt_array')"
        )
        add_entry_to_hdf5(
            out_h5, "ms1_rt_array", ms1_rt, (len(ms1_rt),), "float64", "seconds",
            "RAW-Array of the retention times of MS1 spectra (to plot tic overlay plot in combination with 'ms1_tic_array')"
        )

        # Adding MS1 Map Infos
        add_entry_to_hdf5(
            out_h5, "ms1_map_rt_array", ms1_map_rt, (len(ms1_map_rt),), "float64", "none",
            "Transformed RT-Array of all single MS1 peaks (to plot MS1 heatmaps with 'ms1_map_intens_array' and 'ms1_map_mz_array', filtered by "
            + str(args.filter_threshold) + "% of the highest MS1 peak to remove noise)", compression="gzip"
        )
        add_entry_to_hdf5(
            out_h5, "ms1_map_intens_array", ms1_map_intens, (len(ms1_map_intens),), "float64", "seconds",
            "RAW-Array of the retention times of MS2 spectra (to plot MS1 heatmaps with 'ms1_map_rt_array' and 'ms1_map_mz_array', filtered by "
            + str(args.filter_threshold) + "% of the highest MS1 peak to remove noise)", compression="gzip"
        )
        add_entry_to_hdf5(
            out_h5, "ms1_map_mz_array", ms1_map_mz, (len(ms1_map_mz),), "float64", "m/z",
            "RAW-Array of the TICs of MS2 spectra (to plot MS1 heatmaps with 'ms1_map_intens_array' and 'ms1_map_rt_array', filtered by "
            + str(args.filter_threshold) + "% of the highest MS1 peak to remove noise)", compression="gzip"
        )

        # Adding MS2 HeatMap Infos
        add_entry_to_hdf5(
            out_h5, "ms2_tic_array", ms2_tic, (len(ms2_tic),), "float64", "none",
            "RAW-Array of the TICs of the MS2 spectra (to plot MS2 heatmaps with 'ms2_rt_array' and 'ms2_mz_array')", compression="gzip"
        )
        add_entry_to_hdf5(
            out_h5, "ms2_rt_array", ms2_rt, (len(ms2_rt),), "float64", "seconds",
            "RAW-Array of the retention times of MS2 spectra  (to plot MS2 heatmaps with 'ms2_tic_array' and 'ms2_mz_array')", compression="gzip"
        )
        add_entry_to_hdf5(
            out_h5, "ms2_mz_array", ms2_mz, (len(ms2_mz),), "float64", "m/z",
            "RAW-Array of the m/z of MS2 spectra  (to plot MS2 heatmaps with 'ms2_rt_array' and 'ms2_tic_array')", compression="gzip"
        )
        ###

        ### Calculate the diffenrece in RT between 25%, 50%, 75% and 100% of data (relative) (RT-MS1_Q1-4)
        # Normalize data between 0 and 1
        rt_ms1_q25_temp = np.quantile(ms1_rt, 0.25)
        rt_ms1_q50_temp = np.quantile(ms1_rt, 0.50)
        rt_ms1_q75_temp = np.quantile(ms1_rt, 0.75)

        rt_ms1_q025 = (rt_ms1_q25_temp - 0)  / rt_duration
        rt_ms1_q050 = (rt_ms1_q50_temp - rt_ms1_q25_temp) / rt_duration
        rt_ms1_q075 = (rt_ms1_q75_temp - rt_ms1_q50_temp) / rt_duration
        rt_ms1_q100 = (rt_duration - rt_ms1_q75_temp) / rt_duration

        for val, desc, key in zip (
            [rt_ms1_q025 ,rt_ms1_q050 ,rt_ms1_q075 ,rt_ms1_q100],
            ["first", "second", "third", "last"],
            ["RT_MS1_Q_000-025" ,"RT_MS1_Q_025-050" ,"RT_MS1_Q_050-075" ,"RT_MS1_Q_075-100"]
        ):
            add_entry_to_hdf5(
                out_h5, key, val, (1,), "float64", "none",
                "The time interval when the " + desc + " 25 % of the total number of MS1 spectra are measured"
            )
        ###

        ### Calculate the diffenrece in RT between 25%, 50%, 75% and 100% of data (relative) (RT-MS2_Q1-4)
        rt_ms2_q25_temp = np.quantile(ms2_rt, 0.25)
        rt_ms2_q50_temp = np.quantile(ms2_rt, 0.50)
        rt_ms2_q75_temp = np.quantile(ms2_rt, 0.75)

        rt_ms2_q025 = (rt_ms2_q25_temp - 0)  / rt_duration
        rt_ms2_q050 = (rt_ms2_q50_temp - rt_ms2_q25_temp) / rt_duration
        rt_ms2_q075 = (rt_ms2_q75_temp - rt_ms2_q50_temp) / rt_duration
        rt_ms2_q100 = (rt_duration - rt_ms2_q75_temp) / rt_duration

        for val, desc, key in zip (
            [rt_ms2_q025 ,rt_ms2_q050 ,rt_ms2_q075 ,rt_ms2_q100],
            ["first", "second", "third", "last"],
            ["RT_MS2_Q_000-025" ,"RT_MS2_Q_025-050" ,"RT_MS2_Q_050-075" ,"RT_MS2_Q_075-100"]
        ):
            add_entry_to_hdf5(
                out_h5, key, val, (1,), "float64", "none",
                "The time interval when the " + desc + " 25 % of the total number of MS2 spectra are measured"
            )
        ###

        ### Calculate the diffenrece in RT between 25%, 50%, 75% and 100% of data (relative) (RT-TIC_Q1-4)
        ms1_ms2_rt_tic = sorted(zip(ms1_rt + ms2_rt, ms1_tic + ms2_tic), key=lambda x: x[0] )
        ms1_ms2_rt = [x[0] for x in ms1_ms2_rt_tic]
        ms1_ms2_tic = [x[1] for x in ms1_ms2_rt_tic]
        cummulative_tics_up_to_ms1_ms2 = []
        for i in range(1, len(ms1_ms2_rt)):
            cummulative_tics_up_to_ms1_ms2.append(
                sum(ms1_ms2_tic[0:i]) # Get the sum of tics up to a point
            )

        # Normalize data between 0 and 1
        np_cummulative_tics_up_to_ms1_ms2 = np.array(cummulative_tics_up_to_ms1_ms2)
        np_cummulative_tics_up_to_ms1_ms2_relative = np_cummulative_tics_up_to_ms1_ms2 / cummulative_tics_up_to_ms1_ms2[-1]
        
        # Create Masks to bin data in 25% of TIC
        mask_025 = np_cummulative_tics_up_to_ms1_ms2_relative < 0.25
        mask_050 = np_cummulative_tics_up_to_ms1_ms2_relative < 0.5
        mask_075 = np_cummulative_tics_up_to_ms1_ms2_relative < 0.75

        rt_ms1_ms2_part_025 = (ms1_ms2_rt[np.where(mask_025 == False)[0][0] - 1] - ms1_ms2_rt[0]) / rt_duration # Get RT-Difference from first and last entry in first 25%
        rt_ms1_ms2_part_050 = (ms1_ms2_rt[np.where(mask_050 == False)[0][0] - 1] - ms1_ms2_rt[np.where(mask_025 == False)[0][0]]) / rt_duration # Get RT-Difference from first and last entry in 25% - 50%
        rt_ms1_ms2_part_075 = (ms1_ms2_rt[np.where(mask_075 == False)[0][0] - 1] - ms1_ms2_rt[np.where(mask_050 == False)[0][0]]) / rt_duration # Get RT-Difference from first and last entry in 50% - 75%
        rt_ms1_ms2_part_100 = (ms1_ms2_rt[-1] - ms1_ms2_rt[np.where(mask_075 == False)[0][0]]) / rt_duration # Get RT-Difference from first and last entry in 75% - 100%
        
        for val, desc, key in zip (
            [rt_ms1_ms2_part_025 ,rt_ms1_ms2_part_050 ,rt_ms1_ms2_part_075 ,rt_ms1_ms2_part_100],
            ["first", "second", "third", "last"],
            ["RT_TIC_Q_000-025", "RT_TIC_Q_025-050", "RT_TIC_Q_050-075", "RT_TIC_Q_075-100"]
        ):
            add_entry_to_hdf5(
                out_h5, key, val, (1,), "float64", "none",
                "The time interval when the " + desc + " 25% of TIC accumulates"
            )
        ###

        # Get MS1-Freq-Max (in Hz)
        np_ms1_rt = np.array(ms1_rt)
        ms1_num_of_elements = []
        for i in range(0, len(ms1_rt)):
            ms1_num_of_elements.append(
                sum( # Moving Window of 60 seconds, and sum up the number of elements within.
                    np.logical_and(np_ms1_rt >= ms1_rt[i], np_ms1_rt <= (ms1_rt[i] + 60))
                )
            )
        ms1_freq_max = max(ms1_num_of_elements) / 60

        add_entry_to_hdf5(
            out_h5, "MS1_Freq_Max", ms1_freq_max, (1,), "float64", "Hz",
            "Highest frequency of MS1 spectra in Hertz (iow: the highest number of measured MS1 spectra within a 1 minute window divided by 60)"
        )

        # Get MS2-Freq-Max (in Hz)
        np_ms2_rt = np.array(ms2_rt)
        ms2_num_of_elements = []
        for i in range(0, len(ms2_rt)):
            ms2_num_of_elements.append(
                sum( # Moving Window of 60 seconds, and sum up the number of elements within.
                    np.logical_and(np_ms2_rt >= ms2_rt[i], np_ms2_rt <= (ms2_rt[i] + 60))
                )
            )
        ms2_freq_max = max(ms2_num_of_elements) / 60
        add_entry_to_hdf5(
            out_h5, "MS2_Freq_Max", ms1_freq_max, (1,), "float64", "Hz",
            "Highest frequency of MS2 spectra in Hertz (iow: the highest number of measured MS1 spectra within a 1 minute window divided by 60)"
        )
        ###
        
        # MS1-Density-Q1-4 MS2-Density-Q1-4
        ms1_num_peaks = []
        ms2_num_peaks = []
        for spectrum in exp.getSpectra():
            if spectrum.getMSLevel() == 1:
                ms1_num_peaks.append(len(spectrum.get_peaks()[1]))
            elif spectrum.getMSLevel() == 2:
                ms2_num_peaks.append(len(spectrum.get_peaks()[1]))

        ms1_num_peaks = sorted(ms1_num_peaks)
        ms2_num_peaks = sorted(ms2_num_peaks)

        ms1_density_q25 = np.quantile(ms1_num_peaks, 0.25)
        ms1_density_q50 = np.quantile(ms1_num_peaks, 0.50)
        ms1_density_q75 = np.quantile(ms1_num_peaks, 0.75)

        ms2_density_q25 = np.quantile(ms2_num_peaks, 0.25)
        ms2_density_q50 = np.quantile(ms2_num_peaks, 0.50)
        ms2_density_q75 = np.quantile(ms2_num_peaks, 0.75)

        for val, desc, desc_ms_level, key in zip (
            [ms1_density_q25 ,ms1_density_q50 ,ms1_density_q75 ,ms2_density_q25 ,ms2_density_q50 ,ms2_density_q75],
            ["25%", "50%", "75%", "25%", "50%", "75%"],
            [1,1,1,2,2,2],
            ["MS1_Density_Q1" ,"MS1_Density_Q2" ,"MS1_Density_Q3" ,"MS2_Density_Q1" ,"MS2_Density_Q2" ,"MS2_Density_Q3"]
        ):
            add_entry_to_hdf5(
                out_h5, key, val, (1,), "float64", "none",
                desc + "-quantile of the peak counts per spectrum over all MS" + str(desc_ms_level) + " spectra"
            )
        ###

        # MS1-TIC-Change-Q1-4, see paper for explanation (https://doi.org/10.1021/acs.jproteome.6b00028)
        ms1_tic_change = []
        for i in range(1, len(ms1_rt)):
            ms1_tic_change.append(
                abs(ms1_tic[i-1] - ms1_tic[i])
            )

        ms1_tic_change = sorted(ms1_tic_change)

        ms1_change_max = max(ms1_tic_change)
        ms1_change_q25 = np.quantile(ms1_tic_change, 0.25)
        ms1_change_q50 = np.quantile(ms1_tic_change, 0.50)
        ms1_change_q75 = np.quantile(ms1_tic_change, 0.75)

        ms1_change_q2 = math.log(ms1_change_q50 / ms1_change_q25, 2)
        ms1_change_q3 = math.log(ms1_change_q75 / ms1_change_q50, 2)
        ms1_change_q4 = math.log(ms1_change_max / ms1_change_q75, 2)

        add_entry_to_hdf5(
            out_h5, "MS1-TIC-Change-Q2", ms1_change_q2, (1,), "float64", "none",
            "log2(50%-quantile / 25%-quantile of the TIC-change) "
        )
        add_entry_to_hdf5(
            out_h5, "MS1-TIC-Change-Q3", ms1_change_q3, (1,), "float64", "none",
            "log2(75%-quantile / 50%-quantile of the TIC-change)"
        )
        add_entry_to_hdf5(
            out_h5, "MS1-TIC-Change-Q4", ms1_change_q4, (1,), "float64", "none",
            "log2(100%-quantile / 75%- quantile of the TIC-change)"
        )
        ###

        # MS1-TIC-Q1-4, see paper for explanation (https://doi.org/10.1021/acs.jproteome.6b00028)
        ms1_tic_q25 = np.quantile(ms1_tic, 0.25)
        ms1_tic_q50 = np.quantile(ms1_tic, 0.50)
        ms1_tic_q75 = np.quantile(ms1_tic, 0.75)

        ms1_tic_q2 = math.log(ms1_tic_q50 / ms1_tic_q25, 2)
        ms1_tic_q3 = math.log(ms1_tic_q75 / ms1_tic_q50, 2)
        ms1_tic_q4 = math.log(max(ms1_tic) / ms1_tic_q75, 2)

        add_entry_to_hdf5(
            out_h5, "MS1-TIC-Q2", ms1_tic_q2, (1,), "float64", "none",
            "log2(50%-quantile / 25%-quantile of the actual TIC)"
        )
        add_entry_to_hdf5(
            out_h5, "MS1-TIC-Q3", ms1_tic_q3, (1,), "float64", "none",
            "log2(75%-quantile / 50%-quantile of the actual TIC)"
        )
        add_entry_to_hdf5(
            out_h5, "MS1-TIC-Q4", ms1_tic_q4, (1,), "float64", "none",
            "log2(100%-quantile / 75%-quantile of the actual TIC)"
        )
