# Next-QC-Flow
Transform the Quality Control workflow from Knime into a workflow in Nextflow

Ziel: Anzahl MS1, MS2, PSM, RT-time (in Form einer DB) mit visueller QC-Kontrolle    
Input: Raw files, fasta dbs    
Output: Plots + Tabellen + Identifizierung in einem Webserver mit Excel Download    


## Workpackages

### 1. Database (MySQL Communication)
Postgres vs SQlight

### 2. I/O (raw files, db) 
 Exchange? / Backup Server
Bash (als NF Input)
 
### 3. Pre-Ident Thermo Raw file parser)
Konvertierung, Peak Picker, ...
TRFP, OpenMS (gleicher Algo?) MS1 nicht picken?

### 4. Identifizierung (Mascot)
evtl. durch Comet ersetzten, CLI, MSGF+

### 5. Inferenz (PIA)
PIA , CLI

### 6. Check SpikeIns
TRFP (JSON input unklar?) , read intensity [ISA]

### 7. FeatureFinder (MS1 id. Peptides)
OpenMS

### 8. Visualisierung
The visualization is done in Python with the plotly package for interactive plots.

The data stem from the different tables in the database. As parameters only the names of the raw files to plot have to be given. In case of regular QC, groups can be defined for colouring the PCA plots. If use_groups = False, the PCA plots will be coloured by timestamp, so that effects of the run order can be assessed.

List of plots (regular QC):

- Figure 0: Show table with all QC metrics (TODO)
- Figure 1: Barplot for total number of MS1 and MS2 spectra (absolute numbers)
- Figure 2: Barplot for number of PSMs, peptides, proteins and identified features (absolute numbers)
- Figure 3: TIC Overlay as Lineplot
- Figure 4: Barplot TIC quartiles
- Figure 5: Barplot MS1 TIC quartiles
- Figure 6: Barplot MS2 TIC quartiles
- Figure 7: Barplot of precursor charge states (relative)
- Figure 8: Barplot of PSM charge states (of identified spectra) (relative)
- Figure 9: Barplot of missed cleavages of PSMs (relative)
- Figure 10: PCA on all data (a) (+ table and plot for Loadings (b) (to assess importance of variables for the PCA))
- Figure 11: PCA on raw data (a) (+ table and plot for Loadings (b) (to assess importance of variables for the PCA))


List of plots (ISA QC):

- Figure 0: Show table with important ISA QC metrics (TODO, this corresponds to the run_data table)
- Figure 1: Barplot with number of proteins, protein groups and unfiltered protein groups (absolute numbers)
- Figure 2: Barplot with number of peptides (absolute numbers)
- Figure 3: Barplot with number of PSMs (absolute numbers)

## Docker Image

A docker image containing all needed dependencies is provided in `docker`.

To build this image (with docker-buildx) execute the following in the root directory of this repository:

> docker build -t nextqcflow:latest . -f docker/Dockerfile

An interactive shell can be started with: `docker run -it nextqcflow:latest` and all workflows can be executed via `nextflow` if necessary
