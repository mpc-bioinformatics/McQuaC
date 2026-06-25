# nf-core/macproqc: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory set with `--outdir`.

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes mass spectrometry data through the following steps:

- [Spectra preparation](#spectra-preparation) — decompression and conversion to mzML
- [Raw-data metrics](#raw-data-metrics) — spectrum-level QC metrics from mzML _(migrating)_
- [Peptide identification](#peptide-identification) — database search with Comet _(migrating)_
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

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
