#!/usr/bin/env python

import sys
import os
import argparse
import base64
import zipfile
import pyopenms
from collections import defaultdict

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-featurexml", help="FeatureXML with already annotated identifications")
    parser.add_argument("-hills", help="Hills-file generated from Dinosaur")
    parser.add_argument("-out_csv", help="The Output statistics CSV")
    parser.add_argument("-report_up_to_charge", default=5, help="Upper limit of range to be reported in a csv table for the charge")

    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_setup()

    # Count features
    features = pyopenms.FeatureMap()
    pyopenms.FeatureXMLFile().load(args.featurexml, features)

    total_num_features = 0
    num_features_charge = defaultdict(lambda: 0)
    total_num_ident_features = 0
    num_ident_features_charge = defaultdict(lambda: 0) # Identified by charge state
    for f in features:

        # Get data to be able to count
        charge = f.getCharge()
        idents_of_f = f.getPeptideIdentifications()
        is_identified = len(idents_of_f) > 0

        # Count accordingly
        # Count features with charge
        total_num_features += 1
        num_features_charge[charge] += 1

        # Count identified, if any
        if is_identified:
            total_num_ident_features += 1
            num_ident_features_charge[charge] += 1

    zipfile.ZipFile("featurexml.zip", mode="w", compresslevel=9).write(args.featurexml, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9, arcname=args.featurexml.split(os.sep)[-1])
    with open("featurexml.zip", "rb") as fb:
        feature_str_bs64 = base64.b64encode(fb.read())

    zipfile.ZipFile("hills.zip", mode="w", compresslevel=9).write(args.featurexml, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9, arcname=args.hills.split(os.sep)[-1])
    with open("hills.zip", "rb") as fb:
        hills_str_bs64 = base64.b64encode(fb.read())

    with open(args.out_csv, "w") as csv_out:

        # First create header:
        header = ["total_num_features", "total_num_ident_features"] + \
            ["num_features_charge_" + str(i) for i in range(1, int(args.report_up_to_charge) + 1)] + \
            ["num_ident_features_charge_" + str(i) for i in range(1, int(args.report_up_to_charge) + 1)] + \
            ["feature_data.featureXML.zip", "feature_data.hills.csv.zip"]

        row = [str(total_num_features), str(total_num_ident_features)] + \
            [str(num_features_charge[i]) for i in range(1, int(args.report_up_to_charge) + 1)] + \
            [str(num_ident_features_charge[i]) for i in range(1, int(args.report_up_to_charge) + 1)] + \
            [feature_str_bs64.decode(), hills_str_bs64.decode()]

        csv_out.write(
            ",".join(header) + "\n"
        )

        csv_out.write(
            ",".join(row) + "\n"
        )


