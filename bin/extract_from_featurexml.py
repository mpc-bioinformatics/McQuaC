#!/usr/bin/env python

import os
import argparse
import base64
import zipfile
import pyopenms
from collections import defaultdict
import h5py
import numpy as np

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-featurexml", help="FeatureXML with already annotated identifications")
    parser.add_argument("-out_hdf5", help="The Output statistics HDF5")
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

    with h5py.File(args.out_hdf5, 'w') as out_h5:
        # First create header:
        header = ["total_num_features", "total_num_ident_features"] + \
            ["num_features_charge_" + str(i) for i in range(1, int(args.report_up_to_charge) + 1)] + \
            ["num_ident_features_charge_" + str(i) for i in range(1, int(args.report_up_to_charge) + 1)] 

        # Then create data
        row = [str(total_num_features), str(total_num_ident_features)] + \
            [str(num_features_charge[i]) for i in range(1, int(args.report_up_to_charge) + 1)] + \
            [str(num_ident_features_charge[i]) for i in range(1, int(args.report_up_to_charge) + 1)]

        # And also its description
        description = [
            "Total number of features found in raw file.",
            "Total number of features with an annotated identificaiton.",
            *["Total number of features with charge " + str(i)  for i in range(1, int(args.report_up_to_charge) + 1)],
            *["Total number of identified features with charge " + str(i) for i in range(1, int(args.report_up_to_charge) + 1)]
        ]

        for header, row, desc in zip (header, row, description):
            out_h5.create_dataset(header, (1,), dtype="int32")
            out_h5[header].attrs["unit"] = "none"
            out_h5[header].attrs["Description"] = desc
            out_h5[header].write_direct(np.array(row, dtype="int32"))

        # Save binary blob
        out_h5.create_dataset("feature_data.featureXML.zip", data=feature_str_bs64.decode())
        out_h5["feature_data.featureXML.zip"].attrs["unit"] = "ZIP File"
        out_h5["feature_data.featureXML.zip"].attrs["Description"] = "A featureXML file in a ZIP container containing all features. This file was base64 encoded to be represented as a string"
