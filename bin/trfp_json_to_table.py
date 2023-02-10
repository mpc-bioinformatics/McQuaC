#!/bin/env python

import sys
import csv
import json
import argparse

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-itrfp_json", help="Input-JSON generated by the ThermoRawFileParser XIC")
    parser.add_argument("-iassociations", help="Queries we entered into ThermoRawFileParser (same order)")
    parser.add_argument("-ocsv", help="The Output statistics CSV")

    return parser.parse_args()
if __name__ == "__main__":

    args = argparse_setup()

    # args.itrfp_json = "/home/luxii/Downloads/fsatas/OEI37731stdtest.json"
    # args.iassociations = "/home/luxii/git/Next-QC-Flow/association.txt"
    # args.ocsv = "output_stat.csv"

    json_struct = []
    associations = []
    accessions = dict()

    # Load XICs
    with open(args.itrfp_json, "r") as trfp_in:
        xics = json.load(trfp_in)

    # Load associations
    associations = list()
    with open(args.iassociations, "r") as ass_in:
        for l in ass_in:
            cols = l[:-1].split(",")
            associations.append(cols)

    pass


    # Now create final for this file

    headers_included = []
    final_table = dict()
    headers = [
        "run_file" ,
        "FIXED_MPCSPIKE1_PEP_XXX_RT_XXX",
        "IDENT_MPCSPIKE1_COUNT",
        "IDENT_MPCSPIKE1_DELTA_RT",
        "IDENT_MPCSPIKE1_PEP_XXX_RT_DELTA",
    ]

    for ass, xic in zip(associations, xics["Content"]):
        if ass[0] not in headers_included:
            # If not already included add all information including headers
            final_table[
                "FIXED_" + ass[0] + "_PEP_" + ass[1] + "_MZ_" + ass[2] + "_RT_" + ass[3]
            ] = sum(xic["Intensities"])

            final_table[
                "IDENT_" + ass[0] + "_COUNT"
            ] = 0

            final_table[
                "IDENT_" + ass[0] + "_DELTA_RT"
            ] = None

            final_table[
                "IDENT_" + ass[0] + "_PEP_" + ass[1] + "_MZ_" + ass[2] + "_RT_DELTA"
            ] = None

            headers_included.append(ass[0])
        else:
            # We found an Identification:
            final_table[
                "IDENT_" + ass[0] + "_COUNT"
            ] += 1

            # Retrieve the info from previous associations

            for ass_ass in associations:
                if ass_ass[0] == ass[0]:
                    final_table[
                        "IDENT_" + ass[0] + "_DELTA_RT"
                    ] = (float(ass[3]) - float(ass_ass[3]))  # Subtract found_rt by expecte_rt

                    final_table[
                        "IDENT_" + ass[0] + "_PEP_" + ass_ass[1] + "_MZ_" + ass_ass[2] + "_RT_DELTA"
                    ] = sum(xic["Intensities"])
                    break


    with open(args.ocsv, "w") as final_output:
        writer = csv.DictWriter(final_output, fieldnames=final_table.keys())

        writer.writeheader()
        writer.writerow(final_table)
    pass

