# nf-core/macproqc

**nf-core/macproqc** (**Ma**ss **C**entric **Pro**teomics **Q**uality **C**ontrol) is a bioinformatics pipeline for comprehensive quality control of mass spectrometry-based proteomics experiments. It accepts raw instrument files from Thermo Fisher (`.raw`) and Bruker (`.d`) instruments as well as pre-converted mzML files and produces a rich set of QC metrics, interactive visualisations, and standardised [mzQC](https://hupo-psi.github.io/mzQC/) output.

The pipeline performs the following main steps:

- **Spectra preparation** — raw vendor files are decompressed and converted to the open mzML format using ThermoRawFileParser (Thermo) or tdf2mzml (Bruker).
- **Raw-data QC metrics extraction** _(migrating)_ — MS1/MS2 spectrum counts, TIC quartiles, precursor charge distributions, retention-time coverage and more are extracted from the mzML files using pyOpenMS.
- **Peptide identification** _(migrating)_ — MS2 spectra are searched against a user-supplied FASTA database, without and if necessary with label information.
- **FDR filtering and protein inference** _(migrating)_ — PSM-level results are filtered at 1 % FDR and protein groups are inferred with PIA — Protein Inference Algorithms.
- **QC metrics for peptide features** _(migrating)_ — OpenMS identifies isotope features; found features are mapped to identifications to yield identification rates and QC metrics are extracted.
- **mzQC output** _(migrating)_ — all QC metrics are exported in the standardised [mzQC](https://hupo-psi.github.io/mzQC/) format for interoperability.
- **Visualisation** _(migrating)_ — an interactive report with barplots, TIC overlays, ion maps and PCA plots is generated using Plotly.

## Documentation

- [Usage](usage.md) — how to run the pipeline, samplesheet format, and all parameters.
- [Parameters](parameters.md) — full parameter reference generated from the pipeline schema.
- [Output](output.md) — description of all output files and directories.
- [Contributing](CONTRIBUTING.md) — how to contribute to the pipeline.
