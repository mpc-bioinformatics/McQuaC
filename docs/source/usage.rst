Usage
=====

The workflow should be compatible with running nextflow on Linux, MacOS and Windows (WSL2).



.. _requirements:
Requirements
------------

* Nextflow
* Docker
* `make`

Everything else is maintained within Docker containers.

**TL;DR** Many software projects in academia are abandoned or not updated regularly which makes it
difficult to use new and old software in a single native environment as they might need different
version of the same dependency. E.g. thermorawfileparser is actually linked against Mono 5 on Conda
while fisher_py needs Mono 6 to properly work. This might become worse with more dependencies as
the workflow matures. Docker container provide a robust way to keep software within their own
environment while Nextflow manages starting containers, mounting data etc. in a transparent manner.

The use of Docker containers also increases the compatibility between OSs as everything is running
on Linux (VMs) in the background and enables us to run the workflow on one of Nextflow's many
distributed executors (K8s, Apptainer, ...).   

There is also at least on library which includes a vendor library and the permission to distribute
it only via Docker container.


.. _installation:

Installation
------------

Clone the repository

Docker
^^^^^^

`make docker-imgs` to build/pull the necessary Docker images


Run the workflow
""""""""""""""""

```
nextflow run -profile docker main.nf --main_raw_spectra_folder <FOLDER_WITH_RAWS_FILES --main_fasta_file <FASTA_FILE> --main_comet_params <COMET_CONFIG> --main_outdir <FOLDER_FOR_RESULTS>
```

For non-x64-hardware like Apple Silicon (M1, M2) usage
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Does not work.
**TL;DR** Mono is magically detecting it is running on ARM even in Rosetta emulated Docker containers. When correctly entering a x64-specific file it throws a assertion detection as the magically detected host architecture does not match the execution architecture. As Mono is required for ThermoRawFileParser and the Python module `fisher_py` there is currently no chance of getting rid of it. Use a x64 VM in UTM or VirtualBox.

When we get rid of Mono:   
```
env DOCKER_DEFAULT_PLATFORM=linux/amd64 nextflow run -profile docker main.nf --main_raw_spectra_folder <FOLDER_WITH_RAWS_FILES --main_fasta_file <FASTA_FILE> --main_comet_params <COMET_CONFIG> --main_outdir <FOLDER_FOR_RESULTS>
```
**TL;DR** Some of the used software and containers are only available for x64, therefore Docker writes a warning on stderr that another architecture is used instead of the host architecture, which stops the Nextflow execution. Setting the architecture explicitly using the the env var `DOCKER_DEFAULT_PLATFORM` solves the problems.


Run the workflow
""""""""""""""""

```
nextflow run main.nf --main_raw_spectra_folder <FOLDER_WITH_RAWS_FILES --main_fasta_file <FASTA_FILE> --main_comet_params <COMET_CONFIG> --main_outdir <FOLDER_FOR_RESULTS>
```

For non-x64-hardware like Apple Silicon (M1, M2) use
""""""""""""""""""""""""""""""""""""""""""""""""""""

Does not work, even with Rosetta due to mono. See section Docker

That’s it! You should now be able to execute `main.nf` and all the other workflows within.
