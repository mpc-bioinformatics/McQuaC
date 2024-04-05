# Parameter Overview

Overview of all parameters, broken down for each individual step.

## Quick Links to each section

* [File Conversion](#file-conversion)
* [Get Custom Column From File](#get-custom-column-from-file)
* [Get Features in Raws](#get-features-in-raws)
* [Get mzML Chromatogram and more](#get-mzml-chromatogram-and-more)
* [Identificatin via Comet](#identification-via-comet)
* [Pia](#pia)
* [Retrive SpikeIns](#retrieve-spiekeins)

## Short Overview over all parameters

| **Parameter**              | **Short Description**                     | **Default**                |
| -------------------------- | ----------------------------------------- | -------------------------- |
| `--ccff_input_thermo`      | Input Spectra                             | `$PWD/raws`                |
| `--ccff_outdir`            | Output folder for results                 | `$PWD/results`             |
| `--ccff_header_in_raws`    | Headers to extract in Thermo-files        | `<empty>`                  |
| `--ccff_header_in_d`       | Headers to extract in Bruker-files        | See detailed description   |
| `--ccff_header_in_d_names` | Column-Names for extracted Bruker-headers | See detailed description   |
| `--gf_thermo_raws`               | Input mzMLs (centroided?)                         | `$PWD/raws`                            |
| `--gf_ident_files`               | Input mzTabs, filtered by an FDR                  | `$PWD/raws`                            |
| `--gf_outdir`                    | Output folder for results                         | `$PWD/results`                         |
| `--feature_finder_params`        | Set additional parameters to generate features    | `<empty>`                              |
| `--gf_resolution_featurefinder`  | Set resolution the in input spectra for features  | Defaults to high resolution settings   |
| `--gf_considered_charges_low`    | Set lowest considered charge                      | `1`                                    |
| `--gf_considered_charges_high`   | Set highest considered charge                     | `6`                                    |
| `--additional_dinosaur_settings` | Set additional parameters for Dinosaur            | `<empty>`                              |
| `--gf_num_procs_conversion`      | Set the number of parallel processes used         | `MAX_NUM_PROCS`                        |
| `--gmc_thermo_raws`          | Input mzMLs (centroided?) | `$PWD/raws`     |
| `--gmc_outdir`               | Output folder for results | `$PWD/results`  |
| `--gmc_num_procs_conversion` | NOT USED?                 | `MAX_NUM_PROCS` |
| `--ic_raw_folder`                      | Input mzMLs                                            | `$PWD/MZMLs`                                               |
| `--ic_fasta_file`                      | Input FASTA-database                                   | `proteins.fasta`                                           |
| `--ic_search_parameter_file`           | Input search parameters for Comet                      | [this file](../example_configurations/high-high.comet.params) |
| `--ic_outdir`                          | Output folder for results                              | `$PWD/results`                                             |
| `--ic_tda`                             | Flag to decide if a decoy database should be generated | `1`                                                        |
| `--ic_num_parallel_threads_per_search` | Number of parallel processes per search                | `4`                                                        |
| `--pia_idents`          | Input idXML files                         | `$PWD/identifications`                                |
| `--pia_analysis_file`   | Input parameter file for Pia              | [this file](../example_configurations/pia-analysis.json) |
| `--pia_parameters_file` | Input ??? (TODO currently not used)       | `<empty>`                                             |
| `--pia_outdir`          | Output folder for results                 | `$PWD/results`                                        |
| `--pia_memory`          | JVM-parameter, limiting memory allocation | `8g`                                                  |
| `--spk_raw_spectra`          | Input Spectra                               | `$PWD/raws`                                       |
| `--spk_identification_files` | Input mzTABSs, filtered by FDR              | `$PWD/idents`                                     |
| `--spk_spike_ins`            | Input csv table defining spike-Ins          | [this file](../example_configurations/spike_ins.csv) |
| `--spk_outdir`               | Output folder for results                   | `$PWD/results`                                    |
| `--spk_num_procs_extraction` | Number of parallel extractions done at once | `MAX_NUM_PROCS`                                   |


## File Conversion

### Parameter Overview
| **Parameter**                | **Short Description**                         | **Default**   |
| ---------------------------- | --------------------------------------------- | ------------- |
| `--ctm_input_spectra`        | Input Spectra                                 | `$PWD/raws`     |
| `--ctm_outdir`               | Output folder for results                     | `$PWD/results`  |
| `--ctm_additional_params`    | Set additional conversion parameters for TRFP | `<empty>`       |
| `--ctm_num_procs_conversion` | Set the number of parallel conversions        | `MAX_NUM_PROCS` |


### Detailed Description

### `--ctm_input_spectra`

**Description**: Path to the folder, containing the *.raw or *.d measured spectras.

**Possible Values**: `<Path to folder containing .d or .raw files>`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/raw_and_d_files` or `/my/path/to/a/folder`

**References**:


### `--ctm_outdir` 

**Description**: The Output-Folder, where to save the converted spectra.

**Possible Values**: `<Output folder for the converted spectra>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/converted_files` or `/path/to/an/output/folder`

**References**:


### `--ctm_additional_params` 

**Description**: Set additional parameters for the ThermoRawFileParser. See Reference for possible options. Some settings may interfere with already set parameters in this script.

**Possible Values**: See Reference

**Default**: `<empty>`

**Examples**: `-a`

**References**: A complete of possible parameters can be found [here](https://github.com/compomics/ThermoRawFileParser)

**TODO**: Currently not used in process `convert_mzml`  and `convert_thermo_raw_files`

---

## Get Custom Column From File


### Parameter Overview

| **Parameter**              | **Short Description**                     | **Default**                |
| -------------------------- | ----------------------------------------- | -------------------------- |
| `--ccff_input_thermo`      | Input Spectra                             | `$PWD/raws`                |
| `--ccff_outdir`            | Output folder for results                 | `$PWD/results`             |
| `--ccff_header_in_raws`    | Headers to extract in Thermo-files        | `<empty>`                  |
| `--ccff_header_in_d`       | Headers to extract in Bruker-files        | See detailed description   |
| `--ccff_header_in_d_names` | Column-Names for extracted Bruker-headers | See detailed description   |


### Detailed Description

### `--ccff_input_thermo`

**Description**: Path to the folder, containing the *.raw or *.d measured spectras.

**Possible Values**: `<Path to folder containing .d or .raw files>`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/raw_and_d_files` or `/my/path/to/a/folder`

**References**:


### `--ccff_outdir`

**Description**: The Output-Folder, where to save extracted headers from *.d or *.raw files

**Possible Values**: `<Output folde for the extracted headers>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/extracted_headers` or `/my/path/where/to/save/extracted/headers/folder`

**References**:


### `--ccff_header_in_raws`

**Description**: Headers to be extracted from RAW-files. These need to be added as additional parameters prefixed with `-ehtp`/`-thtp`/`-lhtp` (for extra/tune or log headers) followed by the header name as displayed e.g. in FreeStyle. In the example, two headers are extracted. 

**Possible Values**: See Description and Reference

**Default**: `<empty>` (see `thermo_extract_custom_headers.py`, where basic headers are always extracted if available)

**Examples**: `"-lhtp 'Ambient temp. (°C)' –ehtp 'Ion Injection Time (ms)'"` to extract Ambient temperature and Ion Injection time directly from raws

**References**: See FreeStyle or other programs, which can list headers in RAW-files. As an alternative, you can use `thermo_extract_custom_headers.py`, which prints all available headers.


### `--ccff_header_in_d`

**Description**: Headers to be extracted from .d-files. These need to be added as additional parameters (similar as to RAW files). In Bruker-files, there is no distinction between header types and simply `-htp` followed by the header name to be extracted is sufficient.

**Possible Values**: See Description and Reference

**Default**: `"-htp 'Vacuum_CurrentFore' -htp 'Vacuum_Extra4thGauge' -htp 'Vacuum_CurrentHigh' -htp 'Vacuum_CurrentFunnel' -htp 'Digitizer_CurrentTemp' -htp 'TOF_DeviceTempCurrentValue1' -htp 'TOF_DeviceTempCurrentValue2'"` extracting seven headers if available.

**Examples**: `" -htp 'Digitizer_CurrentTemp' "` for extracting the digitizer current temperature only.

**References**: TODO: Provide resource, which lists all available Headers!


### `--ccff_header_in_d_names`

TODO: To Be Removed. This only sets the column names of the above argument. This will be automaized in future

---

## Get Features in Raws


### Parameter Overview
| **Parameter**                    | **Short Description**                             | **Default**                            |
| -------------------------------- | ------------------------------------------------- | -------------------------------------- |
| `--gf_thermo_raws`               | Input mzMLs (centroided?)                         | `$PWD/raws`                            |
| `--gf_ident_files`               | Input mzTabs, filtered by an FDR                  | `$PWD/raws`                            |
| `--gf_outdir`                    | Output folder for results                         | `$PWD/results`                         |
| `--feature_finder_params`        | Set additional parameters to generate features    | `<empty>`                              |
| `--gf_resolution_featurefinder`  | Set resolution the in input spectra for features  | Defaults to high resolution settings   |
| `--gf_considered_charges_low`    | Set lowest considered charge                      | `1`                                    |
| `--gf_considered_charges_high`   | Set highest considered charge                     | `6`                                    |
| `--additional_dinosaur_settings` | Set additional parameters for Dinosaur            | `<empty>`                              |
| `--gf_num_procs_conversion`      | Set the number of parallel processes used         | `MAX_NUM_PROCS`                        |


### Detailed Description

### `--gf_thermo_raws`

**Description**: A path containing (TODO centroided?) mzML files. Their filenames MUST correspond with the mzTAB filenames. These could be generated by the `File Conversion` workflow.

**Possible Values**: `<Path to folder containing mzML files. >`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/converted_files` or `path/to/folder/containing/converted/files`

**References**:

**TODO**: Parameter is not named correctly


### `--gf_ident_files`

**Description**: A path containing the (by an FDR) filtered mzTAB identification files. These names MUST correspond to the mzML files.

**Possible Values**: `<Path to folder already filtered mzTAB identification files >`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/identifications` or `path/to/identiciation/files/folder`

**References**:


### `--gf_outdir`

**Description**: The Output-Folder, where to save the resulting features and statistics.

**Possible Values**: `<Output folder for the statistics and found features>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/feature_statistics` or `/my/path/where/to/save/feature_statistics`

**References**:


### `--feature_finder_params`

**Description**: Additional parameters which are directly passed to the FeatureFinderMultiplex by OpenMS. See reference for parameters which could be set.

**Possible Values**: See Reference

**Default**: `<empty>`

**Examples**: ` -algorithm:rt_min 5 `

**References**: A detailed list of all available parameters can be found [here](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderMultiplex.html)

**TODO**: Parameter is not named correctly


### `--gf_resolution_featurefinder`

**Description**: 

**Possible Values**: See Reference

**Default**: `"-algorithm:mass_trace:mz_tolerance 0.004 -algorithm:isotopic_pattern:mz_tolerance 0.005"`

**Examples**: `"-algortihm:mass_trace:mz_tolerance 0.02 -algorithm:isotopic_pattern:mz_tolerance 0.04"`  (for low resolution Mass Spectrometers)

**References**: A detailed list of all available parameters can be found [here](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/2.7.0/html/TOPP_FeatureFinderCentroided.html)


### `--gf_considered_charges_low`

**Description**: Charges for the feature finder to use to extract features.

**Possible Values**: `<Positive Number (lower then high)>`

**Default**: `1`

**Examples**: `1` (to also search for singly charged features)

**References**:

**TODO**: Only used in FFCentroided!

### `--gf_considered_charges_high`

**Description**: Charges for the feature finder to use to extract features.

**Possible Values**: `<Positive Number (higher then low>`

**Default**: `6`

**Examples**: `10` (discard all features with a charge higher then 10)

**References**:

**TODO**: Only used in FFCentroided!


### `--additional_dinosaur_settings`

**Description**: Additional parameters, directly passed to the Dinosaur JAR. In the reference a list of possible parameters can be found.

**Possible Values**: See Reference

**Default**: `<empty>`

**Examples**: TODO

**References**: The repository lists some parameters [here](https://github.com/fickludd/dinosaur). In source, parameters are described [here](https://github.com/fickludd/dinosaur/blob/master/src/main/scala/se/lth/immun/DinosaurParams.scala) and advanced parameters in the following [file](https://github.com/fickludd/dinosaur/blob/master/src/main/scala/se/lth/immun/DinosaurAdvParams.scala)


### `--gf_num_procs_conversion`

**Description**: Set the number of processes to be used in the FeatureFinderMultiplex in one execution.

**Possible Values**: `<Nonzero Positive Number>`

**Default**: `MAX_NUM_PROCS`

**Examples**: `16`

**References**:

**TODO**: This process only allows one Fork?


## Get mzML Chromatogram and more

### Parameter Overview
| **Parameter**                | **Short Description**     | **Default**     |
| ---------------------------- | ------------------------- | --------------- |
| `--gmc_thermo_raws`          | Input mzMLs (centroided?) | `$PWD/raws`     |
| `--gmc_outdir`               | Output folder for results | `$PWD/results`  |
| `--gmc_num_procs_conversion` | NOT USED?                 | `MAX_NUM_PROCS` |

### Detailed Description

### `--gmc_thermo_raws`

**Description**: A path to the mzML files. These must be peak picked. These could be generated by the `File Conversion` workflow. 

**Possible Values**: `<Path to folder containing mzML files>`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/converted_files` or `path/to/fileconversion/converted/files`

**References**:


### `--gmc_outdir`

**Description**: The Output-Folder, where to save generated mzML statistics.

**Possible Values**: `<Output folder for generated mzML statistics from this workflow>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/mzml_statistics` or `a/path/which/is/use/to/save/mzml/statistics`

**References**:


### `--gmc_num_procs_conversion`

**Description**: 

**Possible Values**: `<Nonzero Positive Number>`

**Default**: `16`

**Examples**: ``

**References**:

**TODO**: This is currently not used in the worklfow. This should limit the number of parallel executions!

---

## Identification via Comet

### Parameter Overview
| **Parameter**                          | **Short Description**                                  | **Default**                                              |
| -------------------------------------- | ------------------------------------------------------ | -------------------------------------------------------- |
| `--ic_raw_folder`                      | Input mzMLs                                            | `$PWD/MZMLs`                                               |
| `--ic_fasta_file`                      | Input FASTA-database                                   | `proteins.fasta`                                           |
| `--ic_search_parameter_file`           | Input search parameters for Comet                      | [this file](../example_configurations/high-high.comet.params) |
| `--ic_outdir`                          | Output folder for results                              | `$PWD/results`                                             |
| `--ic_tda`                             | Flag to decide if a decoy database should be generated | `1`                                                        |
| `--ic_num_parallel_threads_per_search` | Number of parallel processes per search                | `4`                                                        |

### Detailed Description

### `--ic_raw_folder`

**Description**: The path to peak picked mzML files for the Comet-Search-Engine to do identification

**Possible Values**: `<Path to folder containing mzML files>`

**Default**: `$PWD/MZMLs`

**Examples**: `/mnt/data/converted_files` or `/path/to/mzmls/for/search`

**References**:

**TODO**: This parameter is not named correctly!

### `--ic_fasta_file`

**Description**: The path to a FASTA-file. This file should be formatted in the UniProt-format.

**Possible Values**: `<Path to FASTA-file>`

**Default**: `proteins.fasta`

**Examples**: `/mnt/data/proteins.fasta` or `my_custom.fasta`

**References**:


### `--ic_search_parameter_file`

**Description**: Path to the Comet-Parameter Settings file. This should contain the needed Search Parameters, such as variable and fixed modifications

**Possible Values**: `<Path to Comet-Params file>`

**Default**: `${baseDir}/example_configurations/high-high.comet.params`

**Examples**: `/mnt/data/comet.params` or `some/other/path/to.params`

**References**: See the [comet parameters documentation](https://comet-ms.sourceforge.net/parameters/parameters_202101/)

### `--ic_outdir`

**Description**: The Output-Folder, where to save the identification results from Comet

**Possible Values**: `<Output folder for identified Spectra and Statistics>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/identifications` or some other output path

**References**:


### `--ic_tda`

**Description**: Internal identification parameter for Comet, to know if a decoy database should be generated. If set to 1 a decoy database is generated. If 0, a database containing decoys with the prefix `DECOY_` should be provided. This Parameter overwrites the comet-parameters setting “decoy_search”

**Possible Values**: `0` or `1`

**Default**: `1`

**Examples**: `0`

**References**: The overwritten [parameter](https://comet-ms.sourceforge.net/parameters/parameters_202101/decoy_search.php)


### `--ic_num_parallel_threads_per_search`

**Description**: Number of processes to be used to search one mzML-file. This overwrites the `num_threads` in the comet search parameters file.

**Possible Values**: `<Nonzero Positive Number>`

**Default**: `4`

**Examples**: `8`

**References**:

---

## Pia

### Parameter Overview
| **Parameter**           | **Short Description**                     | **Default**                                           |
| ----------------------- | ----------------------------------------- | ----------------------------------------------------- |
| `--pia_idents`          | Input idXML files                         | `$PWD/identifications`                                |
| `--pia_analysis_file`   | Input parameter file for Pia              | [this file](../example_configurations/pia-analysis.json) |
| `--pia_parameters_file` | Input ??? (TODO currently not used)       | `<empty>`                                             |
| `--pia_outdir`          | Output folder for results                 | `$PWD/results`                                        |
| `--pia_memory`          | JVM-parameter, limiting memory allocation | `8g`                                                  |

### Detailed Description


### `--pia_idents`

**Description**: A path containing the idXML identification files. The identification results should be unfiltered and should contain decoys. Pia filters the identification results.

**Possible Values**: `<Path to the folder containing idXML>`

**Default**: `$PWD/identifications`

**Examples**: `/mnt/data/identifications` or any other path containing identification files

**References**:


### `--pia_analysis_file`

**Description**: A parameter file to configure Pia. IMPORTANT: Pia needs to output mzTAB-files for filtered-PSMs in order for further steps to work properly. See reference for possible parameters. Also, the FDR is set here explicitly!

**Possible Values**: `<Path to parameter file>`

**Default**: `${baseDir}/example_configurations/pia-analysis.json`

**Examples**: `/mnt/data/pia_analysis.json`

**References**: [Project](https://github.com/medbioinf/pia)


### `--pia_parameters_file`

TODO: is this the same as above?


### `--pia_outdir`

**Description**: The Output-Folder, saving the filtered identifications. Here, mzTAB-files are saved. If configuring, make sure to set up PIA to also export mzTAB-files.

**Possible Values**: `<Output folder for filtered identifications / mzTABs)`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/filtered_identifications` or `/the/path/where/to/save/filtered/identifications`

**References**:


### `--pia_memory`

**Description**: A parameter, to let the JVM know, how much memory it can allocate. Depending on the identifications, this parameter may need to be increased

**Possible Values**: `Xg`

**Default**: `8g`

**Examples**: `16g` to give Pia 16G of RAM

**References**:

---

## Retrieve SpiekeIns

### Parameter Overview
| **Parameter**                | **Short Description**                       | **Default**                                       |
| ---------------------------- | ------------------------------------------- | ------------------------------------------------- |
| `--spk_raw_spectra`          | Input Spectra                               | `$PWD/raws`                                       |
| `--spk_identification_files` | Input mzTABSs, filtered by FDR              | `$PWD/idents`                                     |
| `--spk_spike_ins`            | Input csv table defining spike-Ins          | [this file](../example_configurations/spike_ins.csv) |
| `--spk_outdir`               | Output folder for results                   | `$PWD/results`                                    |
| `--spk_num_procs_extraction` | Number of parallel extractions done at once | `MAX_NUM_PROCS`                                   |

### Detailed Description


### `--spk_raw_spectra`

**Description**: Path to the folder, containing the *.raw or *.d measured spectras.

**Possible Values**: `<Path to folder containing .d or .raw files>`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/raw_and_d_files` or any other folder containing raw spectra

**References**:


### `--spk_identification_files`

**Description**: A path containing the (by FDR) filtered mzTAB identification files. These names MUST correspond to the raw spectra files. E.G.: `test1.mzTAB <-> test1.raw` or `buker123.mzTAB <-> bruker123.d`

**Possible Values**: `<Path to folder already filtered mzTAB identification files >`

**Default**: `$PWD/idents`

**Examples**: `/mnt/data/identifications` or any other location

**References**:


### `--spk_spike_ins`

**Description**: A path to a file which should contain: “accession”. “sequence”, “mz” and “RT” as columns. This CSV-file is used to directly extract XICs from the raw-spectra as defined. Additionally, if a identification matches with any of the defined accessions, also xics with the experimentally found mz and RT is extracted.

**Possible Values**: `<Path to CSV-file having accession,sequence,mz,RT as columns>`

**Default**: `${baseDir}/example_configurations/spike_ins.csv`

**Examples**: `/mnt/data/custom_spike_ins.csv`

**References**: Use the example [file](../example_configurations/spike_ins.csv) as a reference


### `--spk_outdir`

**Description**: Output folder of the extracted ion chromatograms (+ some statistical information)

**Possible Values**: `<Output folder for the retrieved SpikeIns statistics)`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/spikein_statistics` or any other output directory

**References**:


### `--spk_num_procs_extraction`

**Description**: This sets the number of parallel XIC-Extractions via TRFP or AlphaTIMS at once. You can limit the number if too many are extracted at once and the machine is running out of memory.

**Possible Values**: `<Nonzero Positive Number>`

**Default**: `MAX_NUM_PROCS`

**Examples**: `16`

**References**:
