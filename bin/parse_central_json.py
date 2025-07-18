#!/usr/bin/env python

import json


def parse_json(json_file_name: str ) -> dict:
    """ Parses the (central) json file and return its content as a dictionary. """
    with open(json_file_name, 'r') as f:
        mcquac_params = json.load(f)

    return mcquac_params


def get_central_param(mcquac_json: dict, param_name: list, default_value=None):
    """ Retrieves a parameter from the central json file's dictionary. The param_name is a list of keys to traverse the dictionary, following"""

    try:
        ret_val = mcquac_json
        for key in param_name:
            ret_val = ret_val[key]

    except KeyError:
        base_peak_tic_up_to = 0
        ret_val = default_value
    
    return ret_val
