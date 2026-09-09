# nf-core/macproqc: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory set with `--outdir`.

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes mass spectrometry data through the following steps:

- [Spectra preparation](#spectra-preparation) — decompression and conversion to mzML
- [Raw-data metrics](#raw-data-metrics) — spectrum-level QC metrics from mzML _(migrating)_
- [Peptide identification](#peptide-identification) — decoy database generation and database search with Comet
- [FDR filtering and protein inference](#fdr-filtering-and-protein-inference) — 1 % FDR filter and PIA protein groups _(migrating)_
- [Feature finding](#feature-finding) — isotope feature detection with IDMapper _(migrating)_
- [mzQC output](#mzqc-output) — standardised QC metric export _(migrating)_
- [Visualization](#visualization) — graphical visualizations of QC metrics using interactive plotly plots
- [Pipeline information](#pipeline-information) — reports from Nextflow itself

---

## Spectra preparation

Raw vendor files are decompressed (if necessary) and converted to the open mzML format. Thermo `.raw` files are converted using ThermoRawFileParser; Bruker `.d` folders are converted using tdf2mzml. Files already in mzML format are passed through after optional decompression.

<details markdown="1">
<summary>Output files</summary>

- `mzmls/`
  - `*.mzML` — converted mzML files. These are intermediate files used by downstream steps and are not published by default

</details>

---

## Raw-data metrics

TODO

---

## Peptide identification

Peptide spectrum matches (PSMs) are produced by searching the converted mzML files against a protein sequence database using the [Comet](https://uwpr.github.io/Comet/) search engine.

### Decoy database

Before the database search, the pipeline generates a target-decoy FASTA by appending decoy sequences (reversed by default) to the input database using the OpenMS `DecoyDatabase` tool. This step can be skipped with `--skip_decoy_generation` if the supplied FASTA already contains decoys.

The decoy database file is an intermediate used downstream for FDR estimation. It is only written to the output directory when `--save_decoy_database true` is set.

<details markdown="1">
<summary>Output files</summary>

- `decoy_database/`
  - `*_decoy.fasta` — combined target-decoy protein sequence database. Only present when `--save_decoy_database true` is set.

</details>

### Comet search configuration

Before the database search, the pipeline generates a Comet parameter file by taking a template (either the built-in default at `assets/default_configs/comet.params` or a user-supplied file via `--comet_config_template`) and overwriting the search settings (mass tolerances, fragment ion parameters, modifications, and output format flags) with the pipeline parameter values.

The adjusted parameter file is an intermediate used by the Comet search step and is not written to the output directory.

### Comet PSMs

Each mzML file is searched against the target-decoy database using Comet and the adjusted parameter file described above. Results are written as mzIdentML (`.mzid`) files.

When `--label_modifications` is set, two searches are run per spectrum file — one unlabelled and one labelled — producing a pair of mzid files per input. The filename encodes the search type:

- `<sample>-unlabelled.mzid` — standard search using `--static_modifications`
- `<sample>-labelled.mzid` — search with `--label_modifications` merged into the static modifications

<details markdown="1">
<summary>Output files</summary>

- `comet/`
  - `*-unlabelled.mzid` — PSM results from the unlabelled search
  - `*-labelled.mzid` — PSM results from the labelled search (only present when `--label_modifications` is set)

</details>

### Visualization

The main output of the visualization subworkflow are interactive plots (as JSON and/or HTML, depending on the settings). Also, a summary output table is generated as csv, tsv or xlsx.

- 00_table_summary.[csv,tsv,xlsx]: summary table containing QC metrics per raw file. The column selection and order is defined by the 'output_column_order' parameter.
- 99_hdf5_feature_table.csv: overview of QC metrics keys and names stored in the hdf5 output files.
- fig01_barplot_MS1_MS2.[plotly.json, html]: Barplot of number of MS1 and MS2 spectra per raw file.
- fig02_barplot_PSMs_peptides_proteins.[plotly.json, html]: Barplot of the number of PSMs, peptides, protein groups and accessions per raw file.
- fig03_barplot_features.[plotly.json, html]: Barplot of the number of features and identified features per raw file.
- fig04_MS1_TIC_overlay.[plotly.json, html]: Lineplot of MS1 Total Ion Chromatogram (TIC) vs retention time for each sample.
- fig05_barplot_TIC_quartiles.[plotly.json, html]: TIC quartiles are the percentage of retention time that is needed to obtain the first, second, third or fourth 25% of the TIC. These quartiles are visualized as stacked barplots.
- fig06_barplot_MS1_TIC_quartiles.[plotly.json, html]: Same as fig05, but only using MS1-level data.
- fig07_barplot_MS2_TIC_quartiles.[plotly.json, html]: Same as fig05, but only using MS2-level data.
- fig08_barplot_precursor_charge.[plotly.json, html]: Stacked barplot of precursor charge fractions (also includes unidentified precursors).
- fig09_barplot_precursor_charge.[plotly.json, html]: Stacked barplot of PSM charge fractions (only identified precursors).
- fig10_barplot_PSM_missedcleavages.[plotly.json, html]: Stacked barplot of missed cleavage fractions (only identified).
- fig11a_PCA_raw.[plotly.json, html]: Principal component analysis plot (PCA-plot). The PCA is calculated on a set of QC metrics that are all directly obtainable from the raw data without a peptide identification step. The metrics are RT_range, nr_MS1, nr_MS2, accumulated_MS1_TIC, accumulated_MS2_TIC, base_peak_intensity_max, total_ion_current_max, MS2_prec_charge_fraction, RT_MS1_quantiles, RT_MS2_quantiles, RT_TIC_quantiles, MS1_freq_max, MS2_freq_max, MS1_density_quantiles, MS2_density_quantiles, MS1_TIC_change_quantiles, MS1_TIC_quantiles.
- fig11b_Loadings_raw.[plotly.json, html]: Loadings plot belonging to the PCA in fig11a. The loadings show which Qc metrics have the highest weight in each direction and can help to interpret groups and outliers in the PCA.
- fig12a_PCA_all.[plotly.json, html]: Principal component analysis plot (PCA-plot). The PCA is calculated on the set of QC metrics used for fig11a plus metrics obtained during peptide identification. The additional metrics are nr_PSMs, nr_peptides, nr_protein_groups, nr_accessions, PSM_charge_fractions, PSM_missed_cleavage_counts, nr_features, nr_ident_features, features_charges, ident_features_charge
- fig12b_Loadings_all.[plotly.json, html]: Loadings plot belonging to the PCA in fig11a.
- fig13_MS1_map: folder containing MS1 maps for each raw file
- fig14_Pump_pressure.[plotly.json, html]: Lineplot of pump pressure vs. retention time for all samples.
- fig15_PSM_error_boxplots.[plotly.json, html]: Boxplots of PSM error (in ppm) for all samples.
- fig16_additional_headers: folder with lineplots over retention time for different metrics recorded by the machine (e.g. temperatures, lock mass correction)
- fig17_BRUKER_calibrant: folder with lineplots over retention time for calibrants used by BRUKER machines

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_<timestamp>.html` — Nextflow execution report with run statistics and per-process resource usage.
  - `execution_timeline_<timestamp>.html` — timeline of all process executions.
  - `execution_trace_<timestamp>.txt` — tab-separated trace file with per-task resource metrics.
  - `pipeline_dag_<timestamp>.html` — directed acyclic graph (DAG) of the pipeline workflow.
  - `nf_core_macproqc_software_versions.yml` — versions of all software used in the pipeline run.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
