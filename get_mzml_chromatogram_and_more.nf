#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.gmc_input_mzmls = "$PWD/mzmls"  // Folder of mzMLs

// Optional Parameters
params.gmc_outdir = "$PWD/results"  // Output-Directory of the mzMLs. Here it is <Input_file>.mzML

/// Tracing show up to 4.7 GB virtual memory for 30000 MS scans
params.get_mzml_chromatogram__get_various_mzml_info = " 7 GB"

// Standalone Workflow
workflow {
    mzmlfiles = Channel.fromPath(params.gmc_mzmls + "/*.mzML")
    get_various_mzml_infos(mzmlfiles)
}


workflow get_various_mzml_infos {
    take:
        mzmlfiles  // mzML where MS1 is peak picked
    main:
        // Retrieve the information from mzml (using pyopenms in a python script)
        retrieve_data_from_mzml(mzmlfiles)
    emit:
        retrieve_data_from_mzml.out
}

process retrieve_data_from_mzml {
    container 'mpc/nextqcflow-python:latest'

    publishDir "${params.gmc_outdir}/", mode:'copy'

    cpus 1
    memory params.get_mzml_chromatogram__get_various_mzml_info

    input:
    path(mzml)

    output:
    path("${mzml.baseName}-mzml_info.csv")

    """
    extract_data_from_mzml.py -mzml ${mzml} -out_csv ${mzml.baseName}-mzml_info.csv
    """
}
