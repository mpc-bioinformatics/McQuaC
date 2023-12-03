#!/usr/bin/env python

import sys
import argparse
import pyopenms
import numpy as np
import math
from collections import defaultdict
import datetime
import time
import csv

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-mzml", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_csv", help="The Output statistics CSV")
    return parser.parse_args()


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
    data_dict = dict()

    # Load MZML
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(args.mzml, exp)

    # Get Date and Time of Measurement (timestamp)
    hh_mm_ss_str = exp.getDateTime().getTime()
    yy_mm_dd_str = exp.getDateTime().getDate() # could also be yy_mm_dd
    try:
        data_dict["timestamp"] = time.mktime(datetime.datetime.strptime(yy_mm_dd_str + "|" + hh_mm_ss_str , "%Y-%m-%d|%H:%M:%S").timetuple())
    except:
        data_dict["timestamp"] = None

    # Get accumulated TICs for ms1 and ms2 (accumulated-MS1_TIC)
    accum_tic_ms1 = get_accumulated_TIC(exp, 1)
    accum_tic_ms2 = get_accumulated_TIC(exp, 2)
    data_dict["accumulated-MS1_TIC"] = accum_tic_ms1
    data_dict["accumulated-MS2_TIC"] = accum_tic_ms2


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

    data_dict["total_num_ms1"] = num_ms1_spectra
    data_dict["total_num_ms2"] = num_ms2_spectra
    data_dict["MS2_PrecZ_Unknown"] = 0 if 0 not in num_ms2_prec_charges else num_ms2_prec_charges[0] / num_ms2_spectra
    for i in range(1, 6):
        data_dict["MS2_PrecZ_" + str(i)] = num_ms2_prec_charges[i] / num_ms2_spectra
    
    data_dict["MS2_PrecZ_more"] = 0
    for key in num_ms2_prec_charges.keys():
        if key > 5:
            data_dict["MS2_PrecZ_more"] += num_ms2_prec_charges[key]
    data_dict["MS2_PrecZ_more"] = data_dict["MS2_PrecZ_more"] / num_ms2_spectra
    
    # Here we have the last spectrum in "spectrum". We get the Retention Time Duration from it
    rt_duration = spectrum.getRT() # In Seconds: How long the MS has acquired spectra 
    data_dict["RT_duration"] = rt_duration

    # Get MS1 and MS2 TICs
    ms1_tic = [] # This is the MS1-TIC-BLOB 
    ms1_rt = []  # This is the MS1-TIC-BLOB
    ms2_tic = [] # This is the MS2-TIC-BLOB
    ms2_rt = []  # This is the MS2-TIC-BLOB
    ms2_mz = []  # This is the MS2-MZ-BLOB (needed or the Ion-Map)
    for spectrum in exp.getSpectra():
        if spectrum.getMSLevel() == 1:
            ms1_rt.append(spectrum.getRT())
            ms1_tic.append(sum(spectrum.get_peaks()[1]))
            precs = spectrum.getPrecursors()
        elif spectrum.getMSLevel() == 2:
            ms2_rt.append(spectrum.getRT())
            ms2_tic.append(sum(spectrum.get_peaks()[1]))
            precs = spectrum.getPrecursors()
            if len(precs) != 1:
                raise Exception("Unexpected number of precursors in a MS2-Spectrum")
            else:
                ms2_mz.append(precs[0].getMZ())
    data_dict["ms1_tic_array"] = ms1_tic
    data_dict["ms1_rt_array"]  = ms1_rt
    data_dict["ms2_tic_array"] = ms2_tic
    data_dict["ms2_rt_array"]  = ms2_rt
    data_dict["ms2_mz_array"]  = ms2_mz

    ### Calculate the diffenrece in RT between 25%, 50%, 75% and 100% of data (relative) (RT-MS1_Q1-4)
    # Normalize data between 0 and 1
    rt_ms1_q25_temp = np.quantile(ms1_rt, 0.25)
    rt_ms1_q50_temp = np.quantile(ms1_rt, 0.50)
    rt_ms1_q75_temp = np.quantile(ms1_rt, 0.75)


    rt_ms1_q025 = (rt_ms1_q25_temp - 0)  / rt_duration
    rt_ms1_q050 = (rt_ms1_q50_temp - rt_ms1_q25_temp) / rt_duration
    rt_ms1_q075 = (rt_ms1_q75_temp - rt_ms1_q50_temp) / rt_duration
    rt_ms1_q100 = (rt_duration - rt_ms1_q75_temp) / rt_duration

    data_dict["RT_MS1_Q_000-025"] = rt_ms1_q025
    data_dict["RT_MS1_Q_025-050"] = rt_ms1_q050
    data_dict["RT_MS1_Q_050-075"] = rt_ms1_q075
    data_dict["RT_MS1_Q_075-100"] = rt_ms1_q100
    ###

    ### Calculate the diffenrece in RT between 25%, 50%, 75% and 100% of data (relative) (RT-MS2_Q1-4)
    rt_ms2_q25_temp = np.quantile(ms2_rt, 0.25)
    rt_ms2_q50_temp = np.quantile(ms2_rt, 0.50)
    rt_ms2_q75_temp = np.quantile(ms2_rt, 0.75)


    rt_ms2_q025 = (rt_ms2_q25_temp - 0)  / rt_duration
    rt_ms2_q050 = (rt_ms2_q50_temp - rt_ms2_q25_temp) / rt_duration
    rt_ms2_q075 = (rt_ms2_q75_temp - rt_ms2_q50_temp) / rt_duration
    rt_ms2_q100 = (rt_duration - rt_ms2_q75_temp) / rt_duration

    data_dict["RT_MS2_Q_000-025"] = rt_ms2_q025
    data_dict["RT_MS2_Q_025-050"] = rt_ms2_q050
    data_dict["RT_MS2_Q_050-075"] = rt_ms2_q075
    data_dict["RT_MS2_Q_075-100"] = rt_ms2_q100
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
    
    data_dict["RT_TIC_Q_000-025"] = rt_ms1_ms2_part_025
    data_dict["RT_TIC_Q_025-050"] = rt_ms1_ms2_part_050
    data_dict["RT_TIC_Q_050-075"] = rt_ms1_ms2_part_075
    data_dict["RT_TIC_Q_075-100"] = rt_ms1_ms2_part_100
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
    data_dict["MS1_Freq_Max"] = ms1_freq_max

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
    data_dict["MS2_Freq_Max"] = ms2_freq_max
        
    
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

    data_dict["MS1_Density_Q1"] = ms1_density_q25
    data_dict["MS1_Density_Q2"] = ms1_density_q50
    data_dict["MS1_Density_Q3"] = ms1_density_q75
    data_dict["MS2_Density_Q1"] = ms2_density_q25
    data_dict["MS2_Density_Q2"] = ms2_density_q50
    data_dict["MS2_Density_Q3"] = ms2_density_q75


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

    data_dict["MS1-TIC-Change-Q2"] = ms1_change_q2
    data_dict["MS1-TIC-Change-Q3"] = ms1_change_q3
    data_dict["MS1-TIC-Change-Q4"] = ms1_change_q4



    # MS1-TIC-Q1-4, see paper for explanation (https://doi.org/10.1021/acs.jproteome.6b00028)
    ms1_tic_q25 = np.quantile(ms1_tic, 0.25)
    ms1_tic_q50 = np.quantile(ms1_tic, 0.50)
    ms1_tic_q75 = np.quantile(ms1_tic, 0.75)

    ms1_tic_q2 = math.log(ms1_tic_q50 / ms1_tic_q25, 2)
    ms1_tic_q3 = math.log(ms1_tic_q75 / ms1_tic_q50, 2)
    ms1_tic_q4 = math.log(max(ms1_tic) / ms1_tic_q75, 2)

    data_dict["MS1-TIC-Q2"] = ms1_tic_q2
    data_dict["MS1-TIC-Q3"] = ms1_tic_q3
    data_dict["MS1-TIC-Q4"] = ms1_tic_q4


    # Get Base_peak_Intensity total_ion_current and up to
    ms1_ms2_basepeaks = []
    ms1_ms2_rt = []
    ms1_ms2_tic = []
    for spectrum in exp.getSpectra():
        ms1_ms2_rt.append(spectrum.getRT())
        peaks = spectrum.get_peaks()[1]
        ms1_ms2_basepeaks.append(sum(peaks) if len(peaks) != 0 else 0)
        ms1_ms2_tic.append(sum(spectrum.get_peaks()[1]))

    base_peak_intensity_max = max(ms1_ms2_basepeaks)
    total_ion_curretn_max = max(ms1_ms2_tic)

    data_dict["Base_Peak_Intensity_Max"] = base_peak_intensity_max
    data_dict["Total_Ion_Current_Max"] = total_ion_curretn_max

    # and up to 105 minutes
    break_on = 105*60
    for i in range(len(ms1_ms2_rt)):
        if break_on < ms1_ms2_rt[i]:
            break
    base_peak_intensity_max_up_to_105m = max(ms1_ms2_basepeaks[:i-1])
    total_ion_curretn_max_up_to_105m = max(ms1_ms2_tic[:i-1])

    data_dict["Base_Peak_Intensity_Max_Up_To_105"] = base_peak_intensity_max_up_to_105m
    data_dict["Total_Ion_Current_Max_Up_To_105"] = total_ion_curretn_max_up_to_105m



    # Write Output
    with open(args.out_csv, "w") as csv_out:

        csv_writer = csv.DictWriter(csv_out, fieldnames=data_dict.keys())

        csv_writer.writeheader()
        csv_writer.writerow(data_dict)
