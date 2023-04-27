#!/bin/env python

"""Extract information for the QC_visualization out of the PIA output files.
TODO 	nrProteins 	nrProteingroups_unfiltered, number-filtered-protein-groups

DONE number-filtered-peptides, psmZ1-5, psm_missed 0-3 ,number-filtered-psms,
"""
import pandas as pd
import csv
import argparse
import zipfile
import base64
import os

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pia_peptides", help="peptides.txt from PIA output")
    parser.add_argument("--pia_proteins", help="Proteins.mzID from PIA output")
    parser.add_argument("-pia_PSMs", help="PSM.mzTab from PIA output")
    parser.add_argument("--output", help="path where output should be saved")

    return parser.parse_args()

def run_pia_extraction():
    args = argparse_setup()
    peptide_count = count_nr_filtered_peptides(args.pia_peptides)
    PSM_counts , charge_counts, miss_counts = read_mzTab(args.PSMs)
    dics = [peptide_count, PSM_counts, charge_counts, miss_counts]

  # geschrieben von Dirk (er muss debuggen falls kaputt geht)
    data = {
        key: [value]
        for d in dics
        for key, value in d.items()
    }
    df = pd.DataFrame(data=data)

    zipfile.ZipFile("pia_extractions.zip", mode="w", compresslevel=9).write(args.pia_peptides,
                                                                       compress_type=zipfile.ZIP_DEFLATED,
                                                                       compresslevel=9,
                                                                       arcname=args.pia_peptides.split(os.sep)[-1])
    zipfile.ZipFile("pia_extractions.zip", mode="w", compresslevel=9).write(args.pia_PSMs,
                                                                       compress_type=zipfile.ZIP_DEFLATED,
                                                                       compresslevel=9,
                                                                       arcname=args.pia_PSMs.split(os.sep)[-1])
    zipfile.ZipFile("pia_extractions.zip", mode="w", compresslevel=9).write(args.pia_proteins,
                                                                       compress_type=zipfile.ZIP_DEFLATED,
                                                                       compresslevel=9,
                                                                       arcname=args.proteins.split(os.sep)[-1])
    with open("pia_extractions.zip", "rb") as pia_b:
        pia_str_bs64 = base64.b64encode(pia_b.read())
        df["pia_output.zip"] = pia_str_bs64
        df.to_csv(args.output)


def count_nr_filtered_peptides(file) -> int:
    #needs the PIA peptides.csv output file and counts PEPTIDE occurences
    with open(file) as peps:
        pepsi_count = 0
        pepsi = csv.reader(peps)
        for row in pepsi:
            if(row[0] == "PEPTIDE"):
                pepsi_count += 1
        peptide_count ={}
        peptide_count["number_filtered_peptides"] = pepsi_count
        return peptide_count

def count_nr_Proteins():
    return int

def read_mzTab(file):
    with open(file) as tabbi:
        PSM_infos = {}

        for lines in tabbi:
            if lines.startswith("PSH"):
                split_lines = lines.split("\t")
                print(split_lines[2], split_lines[14], split_lines[22])
            if lines.startswith("PSM"):
                split_lines = lines.split("\t")
                # PSM_ID [2], charge [14], opt_global_missed_cleavages[22]
                PSM_infos[split_lines[2]] = [split_lines[14], split_lines[22]]


        charge_counts_above5 = 0
        charge_counts_1 = 0
        charge_counts_2 = 0
        charge_counts_3 = 0
        charge_counts_4 = 0
        charge_counts_5 = 0


        miss_count_0 = 0
        miss_count_1 = 0
        miss_count_2 = 0
        miss_count_3 = 0
        miss_count_more = 0
        PSM_counts = len(PSM_infos)
        PSM_count = {}
        PSM_count["number_filtered_pms"] = PSM_counts
        for key, values in PSM_infos.items():
            charge = int(values[0])
            if charge > 5:
                charge_counts_above5 += 1
            if charge == 1:
                charge_counts_1 += 1
            if charge == 2:
                charge_counts_2 += 1
            if charge == 3:
                charge_counts_3 += 1
            if charge == 4:
                charge_counts_4 += 1
            if charge == 5:
                charge_counts_5 += 1

            miss = int(values[1])
            if miss == 0:
                miss_count_0 += 1
            if miss == 1:
                miss_count_1 += 1
            if miss == 2:
                miss_count_2 += 1
            if miss == 3:
                miss_count_3 += 1
            if miss > 3:
                miss_count_more += 1

        charge_counts = {"Z1": charge_counts_1, "Z2": charge_counts_2, "Z3": charge_counts_3, "Z4": charge_counts_4, "Z5": charge_counts_5, "Z_more": charge_counts_above5}
        miss_counts = {"missed_0": miss_count_0, "missed_1": miss_count_1, "missed_2": miss_count_2, "missed_3": miss_count_3, "missed_more": miss_count_more}

        return PSM_count, charge_counts, miss_counts


run_pia_extraction()