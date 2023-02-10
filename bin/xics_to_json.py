#!/bin/env python

import sys
import csv
import json
import argparse

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-icsv", help="Input-CSV-File of the SpikeIns")
    parser.add_argument("-iidents", help="Input-Identification-Files in mzTAB format (here the SpikeIns are specially encoded (see FASTA))")
    parser.add_argument("-ojson", help="The JSON-Output-File to be used with the ThermoRawFileParser XIC")
    parser.add_argument("-oassociation", help="A association file which helps mapping back the the JSON-Results from the ThermoRawFileParser")


    return parser.parse_args()

if __name__ == "__main__":

    args = argparse_setup()

    # args.icsv = "parameter_files/spike_ins.csv"
    # args.iidents = "/home/luxii/Downloads/fsatas/OEI37731stdtest.mzTab"
    # args.ojson = "output.json"
    # args.oassociation = "association.txt"

    json_struct = []
    associations = []
    accessions = dict()
    with open(args.icsv, "r") as in_csv_file:
        # Read CSV and get header indicies
        csv_in = csv.reader(in_csv_file)
        header = next(csv_in)
        ac_idx = header.index("accession")
        sq_idx = header.index("sequence")
        mz_idx = header.index("mz")
        rt_idx = header.index("RT")

        # Generate the two output structures
        for l in csv_in:
            # For the json struct we directly do it in the desired TRFP-Format
            d = dict()
            d["mz"] = float(l[mz_idx])
            d["tolerance"] = 10
            d["tolerance_unit"] = "ppm"
            d["rt_start"] = (float(l[rt_idx]) / 60) - 3 
            d["rt_end"] = (float(l[rt_idx]) / 60) + 3 
            json_struct.append(d)

            # Association to the results is then done per entry from the output-JSON
            associations.append((l[ac_idx], l[sq_idx], l[mz_idx], l[rt_idx]))

            # Have a mapping dictionary for identifications
            accessions[l[ac_idx]] = [l[mz_idx], l[rt_idx], l[sq_idx]]


    with open(args.iidents, "r") as in_ident_file:
        rt_idx = 0
        accession_idx = 0
        for l in in_ident_file:
            # Get Header Line
            if l.startswith("PSH"):
                cols = l.split("\t")
                rt_idx = cols.index("retention_time")
                accession_idx = cols.index("accession")

            # Check if we found an psms which actually fits to the SpikeIns
            if l.startswith("PSM"):
                cols = l.split("\t")
                if cols[accession_idx] in accessions.keys():
                    d = dict()
                    d["mz"] = float(accessions[cols[accession_idx]][0])
                    d["tolerance"] = 10
                    d["tolerance_unit"] = "ppm"
                    d["rt_start"] = (float(cols[rt_idx]) / 60) - 3 
                    d["rt_end"] = (float(cols[rt_idx]) / 60) + 3 
                    json_struct.append(d)

                    associations.append((cols[accession_idx], accessions[cols[accession_idx]][2], accessions[cols[accession_idx]][0], cols[rt_idx]))


    # Write Output JSON-File
    with open(args.ojson, "w") as o_json:
        o_json.write(json.dumps(json_struct, indent=4))

    # Write Output Association-File
    with open(args.oassociation, "w") as output_association:
        for entry in associations:
            output_association.write(
                ",".join(entry) + "\n"
            )
