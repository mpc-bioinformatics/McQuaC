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
plotly, scikit, seaborn,...
