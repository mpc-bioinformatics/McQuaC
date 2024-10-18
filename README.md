# McQuaC

[![Documentation Status](https://readthedocs.org/projects/mcquac/badge/?version=latest)](https://mcquac.readthedocs.io/en/latest/?badge=latest)

## About McQUaC

This workflow does quality control for DDA mass spectrometry proteomics data. The workflow can be used in two different modes (internal standard evaluation = `isa`).

For ISA (`is_isa = TRUE`): Additional steps are executed. Spike-Ins are retrieved from the spectra and also used for evaluation

For no ISA (`is_isa = False`): No additional steps are executed.

Some of the retrieved metrics in this workflow are based on the following publication:

```text
Bittremieux W, Meysman P, Martens L, Valkenborg D, Laukens K. Unsupervised Quality Assessment of Mass Spectrometry Proteomics Experiments by Multivariate Quality Control Metrics. J Proteome Res. 2016 Apr 1;15(4):1300-7. doi: 10.1021/acs.jproteome.6b00028. Epub 2016 Mar 18. PMID: 26974716.
```

The [table](docu/QC_columns_overview.csv) shows all columns which can be retrieved by this quality control workflow, including a brief description for each.

## Getting started

You need to run the workflow with nextflow. For more information see the [documentation](http://mcquac.readthedocs.io/).

## Usage
The workflow should be compatible with Linux, MacOS and Windows (WSL2).

### Requirements

* Nextflow
* Docker
* `make`

Every thing else is maintained within Docker containers.   
TL;DR Many software projects in academia are abandoned or not updated regularly which makes it difficult to use new and old software in a single native environment as they might need different version of the same dependency. E.g. thermorawfileparser is actually linked against Mono 5 on Conda while fisher_py needs Mono 6 to properly work. This might become worse with more dependencies as the workflow matures. Docker container provide a robust way to keep software within their own environment while Nextflow manages starting containers, mounting data etc. in a transparent manner.   
The use of Docker containers also increases the compatibility between OSs as everything is running on Linux (VMs) in the background and enables us to run the workflow on one of Nextflow's many distributed executors (K8s, Apptainer, ...).   
There is also at least on library which includes a vendor library and the permission to distribute it only via Docker container.


### Installation
Clone the repository

#### Docker
`make docker-imgs` to build/pull the necessary Docker images


##### Run the workflow
```
nextflow run -profile docker main.nf --main_raw_spectra_folder <FOLDER_WITH_RAWS_FILES --main_fasta_file <FASTA_FILE> --main_comet_params <COMET_CONFIG> --main_outdir <FOLDER_FOR_RESULTS>
```

##### For non-x64-hardware like Apple Silicon (M1, M2) use
Does not work.
TL;DR Mono is magically detecting it is running on ARM even in Rosetta emulated Docker containers. When correctly entering a x64-specific file it throws a assertion detection as the magically detected host architecture does not match the execution architecture. As Mono is required for ThermoRawFileParser and the Python module `fisher_py` there is currently no chance of getting rid of it. Use a x64 VM in UTM or VirtualBox.

When we get rid of Mono:   
```
env DOCKER_DEFAULT_PLATFORM=linux/amd64 nextflow run -profile docker main.nf --main_raw_spectra_folder <FOLDER_WITH_RAWS_FILES --main_fasta_file <FASTA_FILE> --main_comet_params <COMET_CONFIG> --main_outdir <FOLDER_FOR_RESULTS>
```
TL;DR Some of the used software and containers are only available for x64, therefore Docker writes a warning on stderr that another architecture is used instead of the host architecture, which stops the Nextflow execution. Setting the architecture explicitly using the the env var `DOCKER_DEFAULT_PLATFORM` solves the problems.


##### Run the workflow
```
nextflow run main.nf --main_raw_spectra_folder <FOLDER_WITH_RAWS_FILES --main_fasta_file <FASTA_FILE> --main_comet_params <COMET_CONFIG> --main_outdir <FOLDER_FOR_RESULTS>
```

##### For non-x64-hardware like Apple Silicon (M1, M2) use
Does not work, even with Rosetta due to mono. See section Docker

That’s it! You should now be able to execute `main.nf` and all the other workflows within.

### Arguments

#### Main DDA workflow
| name | default | used in | description |
| --- | --- | --- | --- |
| `--main_raw_spectra_folder` | n/a | `main.nf` | Folder where your raw-files (Thermo) or .d-folders (Bruker) are located |
| `--main_fasta_file` | n/a | `main.nf` | FASTA file containing your targets and optionally spikein if needed |
| `--main_comet_params` | `"${baseDir}/example_configurations/high-high.comet.params"` | `main.nf` | Comet paramter file adjusted to your MS and experiment. See [Comet's documentation](https://uwpr.github.io/Comet/parameters/parameters_202401/) for more information. Use `make comet-params` to generate a new parameters file for the used Comet version. |
| `--spike_ins_table` | `"${baseDir}/example_configurations/spike_ins.csv"` | `main.nf` | CSV with your spike-in peptides information |
| `--main_outdir` | `"$PWD/results"` | `main.nf` | Directory to store the results |
| `--search_spike_ins` | `true` | `main.nf` | Set to false if you do not have any spike-in peptides in your sample |
| `--search_labelled_spikeins` | `true` | `main.nf` | Set to false, if your spike ins are not labelled |

#### Main DIA workflow
Not ready yet

| name | default | used in | description |
| --- | --- | --- | --- |
| `--main_raw_spectra_folder` | `""` | `main_dia.nf` |  |
| `--main_outdir` | `"$PWD/results"` | `main_dia.nf` |  |
| `--main_is_isa` | `true` | `main_dia.nf` |  |

#### Subworkflows
Each subworklow has parameters of its own which are mainly used for scaling. These can be adjusted by simply adding the parameter to the main workflow.

Some of the paramters are not a hard limit but a estimation of how much ressources are used so Nextflow can anticipate how much processes to start. This is especially true for memory usage.


#### Subworkflow: Feature detection
| name | default | used in | description |
| --- | --- | --- | --- |
| `--openms_threads` | `8` | `src/feature_detection.nf` | Max. cores to use per feature finder execution |
| `--min_charge` | `2` | `src/feature_detection.nf` | Minimum charge to consider |
| `--max_charge` | `5` | `src/feature_detection.nf` | Maximum charge to consider |

#### Subworkflow: Comet identification
| name | default | used in | description |
| --- | --- | --- | --- |
| `--identification__comet_threads` | `8` | `src/identification/comet.nf` | Max. cores to use for Comet |
| `--identification__comet_mem` | `"10 GB"` | `src/identification/comet.nf` | Max. memory usage anticipated when running Comet. This is not a hard limit. |

#### Subworkflow: File conversions
| name | default | used in | description |
| --- | --- | --- | --- |
| `--file_conversion__thermo_raw_conversion_mem` | `"10 GB"` | `src/io/raw_file_conversion.nf` | Maximum memory usage anticipated when converting a Thermo rawfile. This is not a hard limit. |
| `--file_conversion__bruker_raw_conversion_mem` | `"5 GB"` | `src/io/raw_file_conversion.nf` | Maximum memory usage anticipated when converting a Bruker .d-folders. This is not a hard limit. |

#### Subworkflow: Inference with PIA
| name | default | used in | description |
| --- | --- | --- | --- |
| `--pia_gb_ram` | `16` | `src/pia.nf` | Maximum memory used by the Java VM and PIA |
| `--pia_threads` | `8` | `src/pia.nf` | Maximum cores used by PIA |

#### Subworkflow: Spike in retrieval
| name | default | used in | description |
| --- | --- | --- | --- |
| `--max_parallel_xic_extractors` | `available cores / 2` | `src/retrieve_spike_ins.nf` |  |

#### Subworkflow: MS metadata extraction
| `--ms_run_metrics__thermo_headers` | n/a | `src/metrics/ms_run_metrics.nf` | Set iuf you want only a specific set of headers extracted |
| `--ms_run_metrics__bruker_headers` | n/a | `src/metrics/ms_run_metrics.nf` | Set iuf you want only a specific set of headers extracted |
| `--ms_run_metrics__thermo_raw_mem` | `"10 GB"` | `src/metrics/ms_run_metrics.nf` | Maximum memory usage anticipated when extracting headers. Not a hard limit. |
| `--ms_run_metrics__bruker_raw_mem` | `"1 GB"` | `src/metrics/ms_run_metrics.nf` | Maximum memory usage anticipated when extracting headers. Not a hard limit. |
| `--ms_run_metrics__mzml_mem` | `"7 GB"` | `src/metrics/ms_run_metrics.nf` | Maximum memory usage anticipated when extracting headers. Not a hard limit. |

## Output

Output is saved in a results folder for each individual file, including all of its step results. The folder `qc_results` within, contains a CSV-table, summarising all generated results. NOTE: This table contains binary blobs and might yield errors, if opening with standard tools, like Excel or LibreOffice. Generated plots (as plotly-html and json)  are also located in this folder.

TODO: The results are saved in a database, which can be later retrieved without executing this workflow again. Please refer to the other provided nextflow-workflows

### Need a GUI?

#### Local
For local use basic GUI is provided. Just install conda and build the    
TODO: Instructions   
activate the environment and run
```
streamlit run simple_main_gui.py
```

#### Want to make it available for your Lab?
Have a look into [MaCWorP](https://github.com/cubimedrub/macworp)! We already added a configuration for integrating McQuaC under `./example_configurations/macworp.json`.

## Detailed General Workflow

- **Data input**: TODO
- **Raw data metrics**: Raw files are converted to mzML, then different metrics are collected and calculated (e.g., number of MS1 and MS2 spectra, TICs, precuror charges).
- **Peptide identification**: The MS2 spectra in the raw files are identified using the Comet search engines using the provided search parameters file and FASTA database.
- **FDR filtering and protein inference**: The results from Comet are further analysed by PIA (filtering for FDR and performing protein inference). PIA creates a mzTAB file, that is further processed.
- **Spike-ins**: Only for ISAs: Identifications and XIC-extraction for a specific list of spike-ins are done on searched spectra files
- **Feature finding**: Using OpenMS FeatureFinder to find features and get metrics like total and identified number of features.
- **TODO: Saving metrics to database**: The output table is saved to the database.
- **Plotting**: In a final step, information is visualized from the generated data and interactive plots are created (like barplots for various metrics, PCA plots or TIC overlay).

![workflow](docu/Workflow_Foto.jpg)

TODO: make the workflow graph nicer

## Details

### 1) Data input

### 2) Raw data metrics

Sub-workflow: `get_mzml_chromatogram_and_more.nf`

The raw files are first converted to mzML using the ThermoRawFileParser with peak picking. The Python script extract_data_from_mzml.py uses the pyopenms package to calculate various metrics like number of MS1 and MS2 spectra, information on retention time and TIC, precursor charge states etc. For details on the different metrics see the [QC_columns_overview.csv](docu/QC_columns_overview.csv) table.

### 3) Peptide identification

Sub-workflows: `identification_via_comet.nf` (TODO and `identification_via_msgfplus.nf`)

The MS2 spectra in the raw files are identified via the Comet (and MSGF+) search engine using the provided FASTA file.

For Comet, a configuration file (`example_configurations/high-high.comet.params`) is needed, which contains the necessary parameters.
By default, peptide mass tolerance is set to 5ppm, fragment tolerance as 20ppm, maximum 2 missed cleavages and Trypsin as the enzyme.
Oxidation(M) is set as a variable modification and Cabamidomethylation (C) as a fixed modification.

**TODO** For MSGF+, the configuration file example_configurations/msgfplus_config.txt is used. Raw files have to be converted to MGF files via the ThermoRawFileParser before identification.
By default, instrument type is set to high-resolution, fragment tolerance as 20ppm, unlimited number of missed cleavages and Trypsin as the enzyme.
Oxidation(M) is set as a variable modification and Cabamidomethylation (C) as a fixed modification.

### 4) FDR filtering and protein inference

Sub-Workflow: `pia.nf`

For FDR estimation and protein inference, the CLI version of PIA (Protein Inference Algorithms) is used. PIA uses an Occam’s Razor principle to combine unique and shared peptides to protein groups. It also applies an FDR filter of 0.01 and contaminants are removed. PIA first has to compile a compilation XML with the input files from the search machines (pia_compilation.xml) which is then used to run the PIA analysis, which results in the identification of PSMs, peptides and proteins.

Input:

The idXML files from the identification steps of the search machines (comet) are used.

Output:

The PIA results are then given as a piaExport-PSMs.mzTab, piaExport-peptides.csv and piaExport-proteins.mzid (though outputs for peptides and proteins might change to mzTab in the future).

Information from the PIA output is extracted in `extract_from_pia_output.py`

From these files, also metrics like the number of identified proteins, proteingroups, PSM charges and missed cleavages are calculated.

### 5) Spike-ins

Sub-Workflow: `retrieve_spike_ins_thermorawfileparser.nf`

For ISA runs only, metrics for a set of predefined spike-in peptides can be retrieved. For this, the spike-in peptides have to be included in the corresponding FASTA file. The list of currently used spike-ins is given in `spike_ins.csv` and can be adapted for other internal standards. We use the ThermoRawFileParser (for .RAW) or AlphaTIMs (for .d) to retrieve the XICs within 10ppm m/z window and a +/- 3 min retention time window of the expected spike-in values.
From the output we aquire the following information for each spike-in:

- `FIXED_MPCSPIKEX_PEP_XXX_RT_XXX` -> Sum of all peaks measured inside the expected spike-in window (10 ppm m/z-tolerance and +/- 3 min retention time window).
- `IDENT_MPCSPIKEX_COUNT` -> spectral count for the given peptide sequence.
- `IDENT_MPCSPIKEX_DELTA_RT` -> difference between expected and measured retention time for the specific spike-in (if peptide was identified).
- `IDENT_MPCSPIKEX_PEP_XXX_RT_DELTA` -> Sum of XIC for the identified PSMs.

### 6) Feature Finding

Sub-workflow. `get_features_in_raws.nf`

The OpenMS FeatureFinderCentroided is used to find the features. The OpenMS IDmapper is used to match the found features with the identifications (mzTabs acquired from PIA). The resulting feature XML file is then used to calculate further metrics like total number of features, charge state of features etc.

### 7) Saving results to database

TODO

### 8) Plotting

The visualization is done in Python with the plotly package for interactive plots. The output csv from the previous step is imported, which may include several rows, one for each raw file.

Input parameters:

| Parameter          | Description                                                                                                                                                                                             |
|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| csv_file           | Path to the CSV file containing the QC results, one row per raw file.                                                                                                                                   |
| group              | List of the experimental group for each sample, comma-separated. Is only used for colouring the PCA plots based on the group. If None, the PCA plots will be coloured by the timestamp of the raw file. |
| output             | Output folder to save the plots.                                                                                                                                                                        |
| fig_show           | Show the figures, only for debugging.                                                                                                                                                                   |
| isa                | True, if it is an ISA run. For ISA only a selected subset of plots are generated.                                                                                                                       |
| tic_overlay_offset | Offset for the TIC plots (TIC of the different raw files are shifted by the offset, default is 0).                                                                                                      |

The plots are saved as `.json` files, which can then be then opened interactively with other tools.

Overview of generated plots:

- Figure 0: Show table with most important QC metrics (**TODO**).
- Figure 1: Barplot for total number of MS1 and MS2 spectra (absolute numbers).
- Figure 2: Barplot for number of PSMs, peptides, protein groups, all filtered (absolute numbers)
- Figure 3: and identified features (absolute numbers).
- Figure 4: Lineplot with the TIC overlay shifted by tic_overlay_offset.
- Figure 5: Barplot of the TIC quartiles.
- Figure 6: Barplot of the MS1 TIC quartiles.
- Figure 7: Barplot of the MS2 TIC quartiles.
- Figure 8: Barplot of the precursor charge states (relative).
- Figure 9: Barplot of PSM charge states (of identified spectra) (relative).
- Figure 10: Barplot of missed cleavages of PSMs (relative).
- Figure 11: PCA on all data (a) (+ plot for Loadings (b) (to assess importance of variables for the PCA)).
- Figure 12: PCA on raw data (a) (+ plot for Loadings (b) (to assess importance of variables for the PCA)).
- Figure 13: Ion map for each raw files (2D density plot with x = RT, y = MZ). These are saved in a separate folder called fig12_ionmaps.
- Figure 14: Time vs. pump pressure (Lineplot)
- Figure 15: Time vs. ion injection time (Lineplot) (still weird)
- Figure 16: Time vs. Lock mass correction (ppm) (Lineplot) (scale is still weird)

In the [column overview](docu/QC_columns_overview.csv) it is indicated which columns are used inside the two PCA plots.

List of plots (ISA QC, TODO!):

- Figure 0: Show table with important ISA QC metrics (TODO, this corresponds to the run_data table)
- Figure 1: Barplot with number of proteins, protein groups and unfiltered protein groups (absolute numbers)
- Figure 2: Barplot with number of peptides (absolute numbers)
- Figure 3: Barplot with number of PSMs (absolute numbers)


## TODO
1. Cleanup all commented code
2. Modularize properly e.g. 
   1. `feature_finding.nf`
      1. `get_features_from_raw()`
   2. `identification.nf`
      1. `comet()`
      2. `msgfplus()`
   3. `inference.nf`
      1. `pia()`
   4. find some name for '...custum_columns...', '...mzml_chromatograms_and_more...', e.g. `sensors.nf`?
3. Create proper module (NF), process (NF) and python script documentation