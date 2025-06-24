#!/usr/bin/env python

import argparse
import json


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-json_in", required=True)
    parser.add_argument("-comet_params", required=True)
    parser.add_argument("-params_out", required=True)
    parser.add_argument("-search_labelled", required=False, default=False, type=bool)
    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    with open(args.json_in, 'r') as f:
        set_params = json.load(f)
    
    if "comet" not in set_params.keys():
        set_params["comet"] = {}

    if "labelled_mods" not in set_params.keys():
        set_params["labelled_mods"] = {}
    
    with open(args.comet_params, 'r') as params_in, open(args.params_out, 'w') as params_out:
        for line in params_in:
            line = line.strip()
            if line.startswith("decoy_search"):                 # disable decoy search
                line = "decoy_search = 0"
            elif line.startswith("output_sqtfile"):             # disable SQT output
                line = "output_sqtfile = 0"
            elif line.startswith("output_txtfile"):             # disable TXT output
                line = "output_txtfile = 0"
            elif line.startswith("output_pepxmlfile"):          # disable pepXML output
                line = "output_pepxmlfile = 0"
            elif line.startswith("output_mzidentmlfile"):       # enable mzIdentML output
                line = "output_mzidentmlfile = 1"
            elif line.startswith("output_percolatorfile"):      # disable percolator output
                line = "output_percolatorfile = 0"
            
            # set the comet params
            for key, value in set_params["comet"].items():
                if line.startswith(key):
                    line = key + " = " + str(value)
            
            # set label parameters (which are set separately, not in comet)
            if args.search_labelled:
                for key, value in set_params["labelled_mods"].items():
                    if line.startswith(key):
                        line = key + " = " + str(value)
            
            params_out.write(line + "\n")
