#!/usr/bin/env python

# std imports
import argparse
from pathlib import Path
import re

PEPTIDE_MASS_TOLERANCE_REGEX = re.compile(r"peptide_mass_tolerance_upper = (?P<tolerance>\d+\.\d+)")
"""Regex to extract the peptide mass tolerance from comet.params
"""

PEPTIDE_MASS_TOLERANCE_UNIT_REGEX = re.compile(r"peptide_mass_units = (?P<unit>\d+)")
"""Regex to extract the peptide mass tolerance unit from comet.params
"""

COMET_FEATURE_FINDER_UNIT_MAP = {
    "0": "DA",
    "2": "ppm"
}
"""Mass units of comet mapped to mass units of the feature finder. 
'2' in comet would be 'mmu' but is not supported by OpenMS's FeatureFinderMultiplex.
"""


def main():
    """Reads mass spectrometer relevant parameters from comet.params and prints them in the format of OpenMS FeatureFinderMultiplex
    """

    # Build and parse CLI
    parser = argparse.ArgumentParser(
        description='Reads mass spectrometer relevant parameters from comet.params and prints them in the format of OpenMS FeatureFinderMultiplex'
    )
    parser.add_argument("-c", "--comet-params", help="The comet params file")
    args = parser.parse_args()

    # Read comet params
    comet_params: str = Path(args.comet_params).read_text(encoding="utf-8")

    # Extract mass tolerance and unit
    precursor_tolerance = PEPTIDE_MASS_TOLERANCE_REGEX.search(comet_params).group("tolerance")
    precursor_tolerance_unit = PEPTIDE_MASS_TOLERANCE_UNIT_REGEX.search(comet_params).group("unit")

    # Map unit to feature finder unit
    feature_finder_unit = COMET_FEATURE_FINDER_UNIT_MAP[precursor_tolerance_unit]

    # Print feature finder params
    print(f"-algorithm:mz_tolerance {precursor_tolerance} -algorithm:mz_unit {feature_finder_unit}")


if __name__ == "__main__":
    main()