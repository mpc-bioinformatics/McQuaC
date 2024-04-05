# Parameter Overview


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

**Possible Values** `<Path to folder containing .d or .raw files>`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/raw_and_d_files` or `/my/path/to/a/folder`

**References**:


### `--ctm_outdir` 

**Description**: The Output-Folder, where to save the converted spectra.

**Possible Values** `<Output folder for the converted spectra>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/converted_files` or `/path/to/an/output/folder`

**References**:


### `--ctm_additional_params` 

**Description**: Set additional parameters for the ThermoRawFileParser. See Reference for possible options. Some settings may interfere with already set parameters in this script.

**Possible Values** See Reference

**Default**: `<empty>`

**Examples**: `-a`

**References**: A complete of possible parameters can be found [here](https://github.com/compomics/ThermoRawFileParser)



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

**Possible Values** `<Path to folder containing .d or .raw files>`

**Default**: `$PWD/raws`

**Examples**: `/mnt/data/raw_and_d_files` or `/my/path/to/a/folder`

**References**:


### `--ccff_outdir`

**Description**: The Output-Folder, where to save extracted headers from *.d or *.raw files

**Possible Values** `<Output folde for the extracted headers>`

**Default**: `$PWD/results`

**Examples**: `/mnt/data/extracted_headers` or `/my/path/where/to/save/extracted/headers/folder`

**References**:


### `--ccff_header_in_raws`

**Description**: Headers to be extracted from RAW-files. These need to be added as additional parameters prefixed with `-ehtp`/`-thtp`/`-lhtp` (for extra/tune or log headers) followed by the header name as displayed e.g. in FreeStyle. In the example, two headers are extracted. 

**Possible Values** See Description and Reference

**Default**: `<empty>` (see `thermo_extract_custom_headers.py`, where basic headers are always extracted if available)

**Examples**: `"-lhtp 'Ambient temp. (°C)' –ehtp 'Ion Injection Time (ms)'"` to extract Ambient temperature and Ion Injection time directly from raws

**References**: See FreeStyle or other programs, which can list headers in RAW-files. As an alternative, you can use `thermo_extract_custom_headers.py`, which prints all available headers.


### `--ccff_header_in_d`

**Description**: Headers to be extracted from .d-files. These need to be added as additional parameters (similar as to RAW files). In Bruker-files, there is no distinction between header types and simply `-htp` followed by the header name to be extracted is sufficient.

**Possible Values** See Description and Reference

**Default**: `"-htp 'Vacuum_CurrentFore' -htp 'Vacuum_Extra4thGauge' -htp 'Vacuum_CurrentHigh' -htp 'Vacuum_CurrentFunnel' -htp 'Digitizer_CurrentTemp' -htp 'TOF_DeviceTempCurrentValue1' -htp 'TOF_DeviceTempCurrentValue2'"` extracting seven headers if available.

**Examples**: `" -htp 'Digitizer_CurrentTemp' "` for extracting the digitizer current temperature only.

**References**: TODO: Provide resource, which lists all available Headers!


### `--ccff_header_in_d_names`

TODO: To Be Removed. This only sets the column names of the above argument. This will be automaized in future


##


### ``

**Description**: 

**Possible Values** ``

**Default**: ``

**Examples**: ``

**References**:






Parameter	Description	Possible Values	Examples	Default	Reference

| Parameter                  | <div style="width:350px">Description</div>                                                                                                                                                                                                                                                  | Possible Values                              | Examples                                                                                                                                | Default                                                                                                                                                                                                                                                                                                                                      | Reference                                                                                                                                                                                                                                                                  |                                   |
|--------------------------  |--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| `--ccff_input_thermo`      | Path to the folder, containing the *.raw or *.d measured spectras.                                                                                                                                                                                           | \<Path to folder containing .d or .raw files> | `/mnt/data/raw_and_d_files`                                                                                                               | `$PWD/raws`                                                                                                                                                                                                                                                                                                                                    |                                                                                                                                                                                                                                                                            |                              |
| `--ccff_outdir`            | The Output-Folder, where to save extracted headers from *.d or *.raw files                                                                                                                                                                                   | \<Output folde for the extracted headers>     | `/mnt/data/extracted_headers`                                                                                                             | `$PWD/results`                                                                                                                                                                                                                                                                                                                                 |                                                                                                                                                                                                                                                                            |                              |
| `--ccff_header_in_raws`    | Headers to be extracted from RAW-files. These need to be added as additional parameters prefixed with -ehtp/thtp/lhtp (for extra/tune or log headers) followed by the header name as displays e.g. in FreeStyle. In the Example, two headers are extracted.  | See Description and Reference                | `-elhtp 'Ambient temp. (°C)' --ehtp 'Ion Injection Time (ms)'` to extract Ambient temperature and Ion Injection time directly from raws | \<empty>                                                                                                                                                                                                                                                                                                                                      | See FreeStyle for all headers. Run this script once, to get a list of all available headers in the *.raw file. This will also display under which category which header falls into.                                                                                        |                                  |
| `--ccff_header_in_d`       | Headers to be extracted from .d-files. These need to be added as additional parameters (similar as to RAW files). In Bruker-files, there is no distinction between header types and simply `-htp` followed by the header name to be extracted is sufficient. | See Description and Reference                | `-htp 'Digitizer_CurrentTemp'` for extracting digitizer                                                                               | `-htp 'Vacuum_CurrentFore' -htp 'Vacuum_Extra4thGauge' -htp 'Vacuum_CurrentHigh' -htp 'Vacuum_CurrentFunnel' -htp 'Digitizer_CurrentTemp' -htp 'TOF_DeviceTempCurrentValue1' -htp 'TOF_DeviceTempCurrentValue2'`                                                                                                                               |                                                                                                                                                                                                                                                                            | No print of all available headers |
| `--ccff_header_in_d_names` |                                                                                                                                                                                                                                                              | See Description and Reference                | `-cn 'BRUKER_Digitizer_CurrentTemp_pickle_zlib'` to name the column in the results table                                              | `-cn 'BRUKER_Vacuum_CurrentFore_pickle_zlib' -cn 'BRUKER_Vacuum_Extra4thGauge_pickle_zlib' -cn 'BRUKER_Vacuum_CurrentHigh_pickle_zlib' -cn 'BRUKER_Vacuum_CurrentFunnel_pickle_zlib' -cn 'BRUKER_Digitizer_CurrentTemp_pickle_zlib' -cn 'BRUKER_TOF_DeviceTempCurrentValue1_pickle_zlib' -cn 'BRUKER_TOF_DeviceTempCurrentValue2_pickle_zlib'` | You can check the bruker files directly to see which headers are available. This can be done by opening the generated spectra files via SQLite. As an alternative, this script can be called, which will report all available headers, which could be used for extraction. | Können wir entfernen              |


## test

![test](test.csv)