#!/usr/bin/env python

import sys
import argparse
import pyopenms
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
    parser.add_argument("-report_up_to_charge", default=5, help="Upper limit of range to be reported in a csv table for the charge")
    return parser.parse_args()


def add_entry_to_hdf5(
    f, qc_acc: str, qc_short_name: str, qc_name: str, qc_description: str, 
    value, value_shape: tuple, value_type: str, 
    unit_accession: str="", unit_name: str=""
    ): 
    """ Adds an entry into the hdf5 file """
    key = "|".join([qc_acc, qc_short_name])  # ACCESSION|SHORT_DESC
    if value_type in ("str", h5py.string_dtype()):
        ds = f.create_dataset(key, shape=value_shape, dtype=h5py.string_dtype(), compression="gzip")
        ds[:] = value
    else:
        f.create_dataset(key, value_shape, dtype=value_type, compression="gzip")
        if not any([x == 0 for x in value_shape]):
            # Check if any dimension is 0. If so, we do not write data in it (zero lengthed result).
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
    args = argparse_setup()

    # # Load MZML
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(args.mzml, exp)

    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:

        # Get Date and Time of Measurement (timestamp)
        hh_mm_ss_str = exp.getDateTime().getTime()
        yy_mm_dd_str = exp.getDateTime().getDate() # could also be yy_mm_dd

        try:
            add_entry_to_hdf5(
                f = out_h5, 
                qc_acc = "LOCAL:10", 
                qc_short_name = "timestamp", 
                qc_name = "UNIX timestamp", 
                qc_description = "Date and time of the measurement as a UNIX timestamp. It is set to -1 if no data is available.", 
                value = time.mktime(datetime.datetime.strptime(yy_mm_dd_str + "|" + hh_mm_ss_str , "%Y-%m-%d|%H:%M:%S").timetuple()), 
                value_shape = (1,), 
                value_type = "float64", 
                unit_accession = "", 
                unit_name = ""
            )
        except:
            add_entry_to_hdf5(
                f = out_h5, qc_acc = "LOCAL:10", 
                qc_short_name = "timestamp", 
                qc_name = "UNIX timestamp", 
                qc_description = "Date and time of the measurement as a UNIX timestamp. It is set to -1 if no data is available.", 
                value = -1, 
                value_shape = (1,), 
                value_type = "int64", 
                unit_accession = "", 
                unit_name = ""
            )

        # Get accumulated TICs for ms1 and ms2 (accumulated-MS1_TIC)
        add_entry_to_hdf5(
            f = out_h5, qc_acc = "MS:4000029", 
            qc_short_name = "accumulated_MS1_TIC", 
            qc_name = "area under TIC in MS1", 
            qc_description = "The area under the total ion current chromatogram (MS:1000235) of all MS1 spectra.", 
            value = get_accumulated_TIC(exp, 1), 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
        )
        add_entry_to_hdf5(
            f = out_h5, qc_acc = "MS:4000030", 
            qc_short_name = "accumulated_MS2_TIC", 
            qc_name = "area under TIC in MS2", 
            qc_description = "The area under the total ion current chromatogram (MS:1000235) of all MS2 spectra.", 
            value = get_accumulated_TIC(exp, 2), 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
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
            f = out_h5, 
            qc_acc = "MS:4000059", 
            qc_short_name = "nr_MS1", 
            qc_name = "number of MS1 spectra", 
            qc_description = "The number of MS1 events in the run.", 
            value = num_ms1_spectra, 
            value_shape = (1,), 
            value_type = "int32", 
            unit_accession = "UO:0000189", 
            unit_name = "count unit"
        )
        
        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000060", 
            qc_short_name = "nr_MS2", 
            qc_name = "number of MS2 spectra", 
            qc_description = "The number of MS2 events in the run.", 
            value = num_ms2_spectra, 
            value_shape = (1,), 
            value_type = "int32", 
            unit_accession = "UO:0000189", 
            unit_name = "count unit"
        )
        
        
        prec_unknown = 0 if 0 not in num_ms2_prec_charges else num_ms2_prec_charges[0] / num_ms2_spectra 
        precz_more = 0
        for key in num_ms2_prec_charges.keys():
            if key > int(args.report_up_to_charge):
                precz_more += num_ms2_prec_charges[key]
        prec_charge_list = [[num_ms2_prec_charges[i] / num_ms2_spectra] for i in range(1, int(args.report_up_to_charge) + 1)]
        ### for now, unknown charge states are reported as 0 (altough it is not mentioned in the ontology)
        prec_charge_list = [[prec_unknown]] + prec_charge_list + [[precz_more / num_ms2_spectra]]
        add_table_in_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000063", 
            qc_short_name = "MS2_prec_charge_fraction", 
            qc_name = "MS2 known precursor charges fractions", 
            qc_description = ("The fraction of MS/MS precursors of the corresponding charge. " 
                              "The fractions [0,1] are given in the 'Fraction' column, corresponding charges in the 'Charge state' column. "
                              "The highest charge state is to be interpreted as that charge state or higher."
                              ), 
            column_names = [str(i) for i in range(0, int(args.report_up_to_charge) + 2)],
            column_data = prec_charge_list, 
            column_types = ["float64" for i in range(0, int(args.report_up_to_charge) + 2)]
        )
        
        
        # Range of the retention time (first and last spectrum)
        RT_first = spectrum_first = exp.getSpectrum(0).getRT()
        RT_last = spectrum.getRT() # Here we have the last spectrum in "spectrum". 
        RT_range = [RT_first, RT_last]
        
        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000070", 
            qc_short_name = "RT_range", 
            qc_name = "retention time acquisition range", 
            qc_description = "Upper and lower limit of retention time at which spectra are recorded.", 
            value = RT_range, 
            value_shape = (2,), 
            value_type = "float64", 
            unit_accession = "UO:0000010", 
            unit_name = "second"
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
            f = out_h5, 
            qc_acc = "LOCAL:01", 
            qc_short_name = "base_peak_intensity_max", 
            qc_name = "base peak intensity max", 
            qc_description = "The maximum base peak (highest peak in spectrum) across all MS1 and MS2 spectra.", 
            value = base_peak_intensity_max, 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
        )
        add_entry_to_hdf5(
            f = out_h5, qc_acc = "LOCAL:03", 
            qc_short_name = "total_ion_current_max", 
            qc_name = "total ion current max", 
            qc_description = "The maximum of all TICs across MS1 and MS2 spectra.", 
            value = total_ion_current_max, 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
        )



        # and up to 105 minutes
        break_on = args.base_peak_tic_up_to*60
        for i in range(len(ms1_ms2_rt)):
            if break_on < ms1_ms2_rt[i]:
                break
        base_peak_intensity_max_up_to_xm = max(ms1_ms2_basepeaks[:i-1])
        total_ion_current_max_up_to_xm = max(ms1_ms2_tic[:i-1])

        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "LOCAL:02", 
            qc_short_name = "base_peak_intensity_max_up_to_" + str(args.base_peak_tic_up_to), 
            qc_name = "base peak intensity max", 
            qc_description = "The maximum base peak (highest peak in spectrum) across all MS1 and MS2 spectra up to " + str(args.base_peak_tic_up_to) + " minutes.", 
            value = base_peak_intensity_max, 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
        )
        add_entry_to_hdf5(
            f = out_h5, qc_acc = "LOCAL:04", 
            qc_short_name = "total_ion_current_max_up_to_" + str(args.base_peak_tic_up_to), 
            qc_name = "total ion current max", 
            qc_description = "The maximum of all TICs across MS1 and MS2 spectra up to " + str(args.base_peak_tic_up_to) + " minutes.", 
            value = total_ion_current_max, 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
        )


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

        add_table_in_hdf5(
            f = out_h5, 
            qc_acc = "LOCAL:25", 
            qc_short_name = "MS1_TIC", 
            qc_name = "MS1 total ion chromatogram", 
            qc_description = ("MS1 TIC over the retention time. " 
                              "This table can be used to plot the MS1 TIC over time. "
                              ), 
            column_names = ["retention_time", "TIC"],
            column_data = [ms1_rt, ms1_tic], 
            column_types = ["float64", "float64"]
        )
       
        add_table_in_hdf5(
            f = out_h5, 
            qc_acc = "LOCAL:26", 
            qc_short_name = "MS1_map", 
            qc_name = "MS1 ion map data", 
            qc_description = ("MS1 intensity, m/z and retention time. " 
                              "This table can be used to plot the MS1 ion map. "
                              ), 
            column_names = ["retention_time", "mz", "intensity"],
            column_data = [ms1_map_rt, ms1_map_mz, ms1_map_intens], 
            column_types = ["float64", "float64", "float64"]
        )
        
        add_table_in_hdf5(
            f = out_h5, 
            qc_acc = "LOCAL:27", 
            qc_short_name = "MS2_TIC", 
            qc_name = "MS2 total ion chromatogram", 
            qc_description = ("MS2 TIC over the retention time. " 
                              "This table can be used to plot the MS1 TIC over time. "
                              ), 
            column_names = ["retention_time", "TIC"],
            column_data = [ms2_rt, ms2_tic], 
            column_types = ["float64", "float64"]
        )
        

        ### Calculate the difference in RT between 25%, 50%, 75% and 100% of data (relative) (RT-MS1_Q1-4)
        # Normalize data between 0 and 1
        rt_ms1_q25_temp = np.quantile(ms1_rt, 0.25)
        rt_ms1_q50_temp = np.quantile(ms1_rt, 0.50)
        rt_ms1_q75_temp = np.quantile(ms1_rt, 0.75)


        rt_duration = RT_last - RT_first
        rt_ms1_q025 = (rt_ms1_q25_temp - 0)  / rt_duration
        rt_ms1_q050 = (rt_ms1_q50_temp - rt_ms1_q25_temp) / rt_duration
        rt_ms1_q075 = (rt_ms1_q75_temp - rt_ms1_q50_temp) / rt_duration
        rt_ms1_q100 = (rt_duration - rt_ms1_q75_temp) / rt_duration

        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000055", 
            qc_short_name = "RT_MS1_quartiles", 
            qc_name = "MS1 quarter RT fraction", 
            qc_description = "The interval used for acquisition of the first, second, third, and fourth quarter of all MS1 events divided by retention time duration.", 
            value = (rt_ms1_q025,rt_ms1_q050,rt_ms1_q075,rt_ms1_q100), 
            value_shape = (4,), 
            value_type = "float64", 
            unit_accession = "UO:0000191", 
            unit_name = "fraction"
        )
        ###

        ### Calculate the difference in RT between 25%, 50%, 75% and 100% of data (relative) (RT-MS2_Q1-4)
        rt_ms2_q25_temp = np.quantile(ms2_rt, 0.25)
        rt_ms2_q50_temp = np.quantile(ms2_rt, 0.50)
        rt_ms2_q75_temp = np.quantile(ms2_rt, 0.75)

        rt_ms2_q025 = (rt_ms2_q25_temp - 0)  / rt_duration
        rt_ms2_q050 = (rt_ms2_q50_temp - rt_ms2_q25_temp) / rt_duration
        rt_ms2_q075 = (rt_ms2_q75_temp - rt_ms2_q50_temp) / rt_duration
        rt_ms2_q100 = (rt_duration - rt_ms2_q75_temp) / rt_duration

        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000056", 
            qc_short_name = "RT_MS2_quartiles", 
            qc_name = "MS2 quarter RT fraction", 
            qc_description = "The interval used for acquisition of the first, second, third, and fourth quarter of all MS2 events divided by retention time duration.", 
            value = (rt_ms2_q025,rt_ms2_q050,rt_ms2_q075,rt_ms2_q100), 
            value_shape = (4,), 
            value_type = "float64", 
            unit_accession = "UO:0000191", 
            unit_name = "fraction"
        )
        
            
        ###
        
        ### Calculate the difference in RT between 25%, 50%, 75% and 100% of data (relative) (RT-TIC_Q1-4)
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
        
        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000054", 
            qc_short_name = "RT_TIC_quartiles", 
            qc_name = "TIC quarters RT fraction", 
            qc_description = "The interval when the respective quarter of the TIC accumulates divided by retention time duration.", 
            value = (rt_ms1_ms2_part_025,rt_ms1_ms2_part_050,rt_ms1_ms2_part_075,rt_ms1_ms2_part_100), 
            value_shape = (4,), 
            value_type = "float64", 
            unit_accession = "UO:0000191", 
            unit_name = "fraction"
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
            f = out_h5, 
            qc_acc = "MS:4000065", 
            qc_short_name = "MS1_freq_max", 
            qc_name = "fastest frequency for MS level 1 collection", 
            qc_description = "fastest frequency for MS level 1 collection.", 
            value = ms1_freq_max, 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "UO:0000106", 
            unit_name = "hertz"
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
            f = out_h5, 
            qc_acc = "MS:4000066", 
            qc_short_name = "MS2_freq_max", 
            qc_name = "fastest frequency for MS level 2 collection", 
            qc_description = "fastest frequency for MS level 2 collection.", 
            value = ms2_freq_max, 
            value_shape = (1,), 
            value_type = "float64", 
            unit_accession = "UO:0000106", 
            unit_name = "hertz"
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

        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000061", 
            qc_short_name = "MS1_density_quartiles", 
            qc_name = "MS1 density quantiles", 
            qc_description = ("The first to n-th quantile of MS1 peak density (scan peak counts). " 
                              "A value triplet represents the original QuaMeter metrics, the quartiles of MS1 density. " 
                              "The number of values in the tuple implies the quantile mode."), 
            value = (ms1_density_q25,ms1_density_q50,ms1_density_q75), 
            value_shape = (3,), 
            value_type = "int32", 
            unit_accession = "UO:0000189", 
            unit_name = "count unit"
        )
        
        add_entry_to_hdf5(
            f = out_h5, 
            qc_acc = "MS:4000062", 
            qc_short_name = "MS2_density_quartiles", 
            qc_name = "MS2 density quantiles", 
            qc_description = ("The first to n-th quantile of MS2 peak density (scan peak counts). " 
                              "A value triplet represents the original QuaMeter metrics, the quartiles of MS2 density. " 
                              "The number of values in the tuple implies the quantile mode."), 
            value = (ms2_density_q25,ms2_density_q50,ms2_density_q75), 
            value_shape = (3,), 
            value_type = "int32", 
            unit_accession = "UO:0000189", 
            unit_name = "count unit"
        )


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
            f = out_h5, 
            qc_acc = "MS:4000057", 
            qc_short_name = "MS1_TIC_change_quartiles", 
            qc_name = "MS1 TIC-change quartile ratios", 
            qc_description = ("The log ratios of successive TIC-change quartiles. " 
                              "The TIC changes are the list of MS1 total ion current (TIC) value changes from one to the next scan, produced when each MS1 TIC is subtracted from the preceding MS1 TIC. "
                              "The metric's value triplet represents the log ratio of the TIC-change Q2 to Q1, Q3 to Q2, TIC-change-max to Q3."), 
            value = (ms1_change_q2,ms1_change_q3,ms1_change_q4), 
            value_shape = (3,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
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
            f = out_h5, 
            qc_acc = "MS:4000057", 
            qc_short_name = "MS1_TIC_quartiles", 
            qc_name = "MS1 TIC quartile ratios", 
            qc_description = ("The log ratios of successive TIC quartiles. " 
                              "The metric's value triplet represents the log ratios of TIC-Q2 to TIC-Q1, TIC-Q3 to TIC-Q2, TIC-max to TIC-Q3."), 
            value = (ms1_tic_q2,ms1_tic_q3,ms1_tic_q4), 
            value_shape = (3,), 
            value_type = "float64", 
            unit_accession = "", 
            unit_name = ""
        )

