#!/usr/bin/env python

import argparse
import h5py
import numpy as np
import os

def check_if_file_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    if os.path.isfile(s):
        return s
    else:
        raise Exception("File '{}' does not exists".format(s))


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hdf_out_name", help="The output of the combined HDF5 files")
    parser.add_argument("-put_under_subdataset", default=False, action="store_true", help="Flag if data should be added under a subdataset")
    parser.add_argument(
        "files", type=check_if_file_exists, nargs="+",
        help="All files which should be combined"
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_setup()

    if not args.put_under_subdataset:
        # Write everything under "/"
        with h5py.File(args.hdf_out_name, 'w') as out_h5:
            for hdf in args.files:
                with h5py.File(hdf, 'r') as in_h5:
                    for obj in in_h5:
                        in_h5.copy(obj, out_h5)
    else: 
        # Write everything under "/filename/"
        with h5py.File(args.hdf_out_name, 'w') as out_h5:
            for hdf in args.files:
                filename = os.path.basename(hdf).split(".", 1)[0]
                out_h5.create_group("/" + filename)
                with h5py.File(hdf, 'r') as in_h5:
                    for obj in in_h5:
                        in_h5.copy(obj, out_h5["/" + filename])
