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
- [Visualisation](#visualisation) — interactive Plotly report _(migrating)_
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
