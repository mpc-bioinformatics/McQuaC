#!/usr/bin/env python

import sys
import csv
import json
import argparse
import zlib
import pickle
import base64
import pandas as pd

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-itrfp_json", help="Input-JSON generated by the ThermoRawFileParser XIC")
    parser.add_argument("-iassociations", help="Queries we entered into ThermoRawFileParser (same order)")
    parser.add_argument("-ocsv", help="The Output statistics CSV")

    return parser.parse_args()
if __name__ == "__main__":

    args = argparse_setup()

    # args.itrfp_json = "/home/luxii/git/Luxii/Next-QC-Flow/work/6d/7e9463f039a556b16eb8be989e4c89/EX05999std.json"
    # args.iassociations = "/home/luxii/git/Luxii/Next-QC-Flow/work/6d/7e9463f039a556b16eb8be989e4c89/association.txt"
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

    # Example header for MPCSPIKE1 (old-isa) (Replace IDENTIFIER/XXX with specific peptide information)
    # headers = [ 
    #     "run_file" ,
    #     "SPIKE_IDENTIFIER_MZ_XXX_RT_XXX_Maximum_Intensity"
    #     "SPIKE_IDENTIFIER_MZ_XXX_RT_XXX_RT_at_Maximum_Intensity"
    #     "SPIKE_IDENTIFIER_MZ_XXX_RT_XXX_FoundPSMs",
    #     "SPIKE_IDENTIFIER_MZ_XXX_RT_XXX_Delta_to_expected_RT",
    # ]

    for ass, xic in zip(associations, xics["Content"]):

        # Map to original SpikeIn (one of the first entries):
        ass = associations[[x[0] for x in associations].index(ass[0])]

        # Get maximum intensity, its RT as well as the delta from the extracted xic
        max_intens = max(xic["Intensities"]) if xic["Intensities"] else None
        rt_at_max_intens = 60*xic["RetentionTimes"][xic["Intensities"].index(max_intens)] if xic["Intensities"] else None
        rt_delta = float(ass[3]) - rt_at_max_intens  #  "Expected retention time" - "Actual Retention Time of SpikeIn (highest peak)""

        # Append inormation to the table
        final_table[
            "SPIKE_" + ass[0] + "_MZ_" + ass[2] + "_RT_" + ass[3] + "_Maximum_Intensity"
        ] = max_intens

        final_table[
            "SPIKE_" + ass[0] + "_MZ_" + ass[2] + "_RT_" + ass[3] + "_RT_at_Maximum_Intensity"
        ] = rt_at_max_intens

        final_table[
            "SPIKE_" + ass[0] + "_MZ_" + ass[2] + "_RT_" + ass[3] + "_Delta_to_expected_RT"
        ] = rt_delta

        if ass[0] not in headers_included: 
            # CASE: Not included in the final table, therefore we have a Spike-In extracted from the initial table

            final_table[
                "SPIKE_" + ass[0] + "_MZ_" + ass[2] + "_RT_" + ass[3] + "_Found_PSMs"
            ] = 0
            
            headers_included.append(ass[0])
        else:
            # CASRA: We already inlcluded this entry. Therefore, this information must come from a PSM!
            final_table[
                "SPIKE_" + ass[0] + "_MZ_" + ass[2] + "_RT_" + ass[3] + "_Found_PSMs"
            ] +=1  # We count up, because we have two entries of it --> Therefore 1 PSM (for three enries, 2 PSMs, etc...)

    # Pickle the final table and create a dataframe with only one column and one line
    final_table_pickled = base64.b64encode(zlib.compress(pickle.dumps(final_table), level=9)).decode("utf-8")
    final_table_pickled = {"SPIKEINS_____pickle_zlib": final_table_pickled}
    final_table_pickled = pd.DataFrame(final_table_pickled, index=[0])

    # Save to csv
    final_table_pickled.to_csv(args.ocsv, index = False)
