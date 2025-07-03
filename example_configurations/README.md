# Files in this folder

## mcquac_params.json

Central configuration file for McQuaC.

This json dictionary holds central settings for the pipeline, which will be parsed for the spectrum identification as well as for the quantification.

The following parameters are supported:


| comet | the value is a dict, all key-value pairs of the dict are passed to the comet parameters in the `comet.params`, like `"fragment_bin_tol": 0.02`. Some of these parameters (like `peptide_mass_tolerance_upper`) are also used for the feature detection. |
|-----|------|
|labelled_mods | value is a dict, interpreted as the (fixed) modifications (comet.params style) which are added to the identification for the "labelled search", if this is activated. |

## spike-ins.csv

This is a custom format for this QC-Workflow. Which gives us the information which PSMs in the identification have/can be found to also extract the XICs from.

Column-Explanation:

`accession`: The actual accession how we could find this entry in an PSM. *NOTE*: We also have to update the FASTA-file to this peptide with the chosen accession.

`sequence`: The sequence which is should be also present in the FASTA. *NOTE*: "old-isa" does not have a sequence and ist not added into the FASTA, hence we can never find identifications for it!

`mz`: The mz, which will be used for the XIC extraction for this entry. We use 10 ppm as a tolerance here

`RT`: The retention time, where we choose the range: (RT / 60) +-3 for XIC-extraction. *Note*: In case of a found PSM, we also extract the XIC of the found PSMs RT.  This will be then annotated in a delta-shift, calculated as "found_psms_RT - original_RT" 