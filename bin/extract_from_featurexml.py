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
