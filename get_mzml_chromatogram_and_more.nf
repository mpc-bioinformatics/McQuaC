#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.gmc_input_mzmls = "$PWD/mzmls"  // Folder of mzMLs

// Optional Parameters
params.gmc_outdir = "$PWD/results"  // Output-Directory of the mzMLs. Here it is <Input_file>.mzML
params.gmc_num_forks = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

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

    maxForks params.gmc_num_forks

    input:
    path(mzml)

    output:
    path("${mzml.baseName}-mzml_info.csv")

    """
    extract_data_from_mzml.py -mzml ${mzml} -out_csv ${mzml.baseName}-mzml_info.csv
    """
}
