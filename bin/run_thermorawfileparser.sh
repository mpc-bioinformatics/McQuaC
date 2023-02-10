#!/bin/bash
mono $(dirname "$0")/ThermoRawFileParser_v1.4.0/ThermoRawFileParser.exe ${@:1}
