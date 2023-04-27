#!/bin/env python

import sys
import os
import argparse
import time
import csv
csv.field_size_limit(sys.maxsize)
from collections import defaultdict


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_csvs", help="List of input_csvs concatenated via ','")
    parser.add_argument("-out_csv", help="The Output statistics CSV")
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    args.input_csvs = "/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38563std_with_idents_funny_outpautpasudifaj.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38536std_mzml_info.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38536std_features.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38536std_spikeins.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38563std_with_idents_funny_outpautpasudifaj.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38560std_features.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38560std_mzml_info.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38560std_spikeins.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38563std_features.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38563std_mzml_info.csv,/home/luxii/git/Next-QC-Flow/results_OLD/QEXI38563std_spikeins.csv"
    args.out_csv = "test_output.csv" 



    # Get timestamp of the actual identification
    ts = time.time()
    

    data_dict = defaultdict(lambda: dict())

    # Gather all the resources
    for input_csv in args.input_csvs.split(","):

        # Get filename
        base_filename = input_csv.split(os.sep)[-1].split("_____", 1)[0]
        if base_filename + "_" + str(ts) not in data_dict:
            data_dict[base_filename + "_" + str(ts)] = defaultdict(lambda: "")
            data_dict[base_filename + "_" + str(ts)]["filename"] = base_filename

        # Get data (and headers)
        with open(input_csv, "r") as in_csv:
            csv_reader = csv.reader(in_csv)
            header = next(csv_reader)
            data = next(csv_reader)

            # And save in a large dictionary
            for h, d in zip(header, data):
                data_dict[base_filename + "_" + str(ts)][h] = d

            

    # Get all headers and include file_and_analysis_timestamp
    all_keys = [list(x.keys()) for x in data_dict.values()]
    final_rows = [["file_and_analysis_timestamp"] + sorted(list(set([y for x in all_keys for y in x])))]


    for key, val in data_dict.items():
        row = [key]
        for header in final_rows[0][1:]:
            row.append(val[header])
        final_rows.append(row)



    with open(args.out_csv, "w") as csv_out:
        csv_writer = csv.writer(csv_out)
        for row in final_rows:
            csv_writer.writerow(
                row
            )


