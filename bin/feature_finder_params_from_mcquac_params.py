#!/usr/bin/env python

import argparse
import json

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-json_in", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    with open(args.json_in, 'r') as f:
        set_params = json.load(f)

    # Extract mass tolerance and unit
    precursor_tolerance = set_params["comet"]["peptide_mass_tolerance_upper"]
    
    if set_params["comet"]["peptide_mass_units"] == 0:
        feature_finder_unit = "DA"
    elif set_params["comet"]["peptide_mass_units"] == 1:
        precursor_tolerance = float(precursor_tolerance) / 1000.0
        feature_finder_unit = "DA"
    elif set_params["comet"]["peptide_mass_units"] == 2:
        feature_finder_unit = "ppm"

    # Print feature finder params
    print(f"-algorithm:mz_tolerance {precursor_tolerance} -algorithm:mz_unit {feature_finder_unit}")
