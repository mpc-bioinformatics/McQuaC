#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ccff_input_spectra = "$PWD/raws"  // Folder of Thermo-RAW-files

// Optional Parameters
params.ccff_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ccff_header_in_raws = ""  // Headers which we want to extract from thermo

params.ccff_header_in_d = "" // Extract headers from Bruker raw measurements-

// Standalone Workflow
workflow {
    raw_spectra = Channel.fromPath(params.ccff_input_spectra + "/*.raw")
    get_custom_headers(raw_spectra)
}


workflow get_custom_headers {
    take:
        thermo_raw_files
        bruker_raw_files
    main:
        // Extract Information directly from the RAW-file using fisher-py
        thermo_headers = retrieve_custom_headers_from_thermo_raw_files(thermo_raw_files)
        bruker_headers = retrieve_custom_headers_from_bruker_raw_files(bruker_raw_files)
        headers = thermo_headers.concat(bruker_headers)
    emit:
        headers
}

process retrieve_custom_headers_from_thermo_raw_files {
    container 'mpc/nextqcflow-python:latest'
    errorStrategy 'ignore'

    stageInMode 'copy'

    input:
    path raw

    output:
    path "${raw.baseName}-custom_headers.csv"

    """
    # Pythonnet sometimes fails to exit and throws a mono error
    thermo_extract_custom_headers.py -raw ${raw} ${params.ccff_header_in_raws} -out_csv ${raw.baseName}-custom_headers.csv || true

    # Fail Check if no content was written
    if ! [ -s "${raw.baseName}-custom_headers.csv" ];then
        rm ${raw.baseName}-custom_headers.csv
    fi
    """
}

process retrieve_custom_headers_from_bruker_raw_files {
    container 'mpc/nextqcflow-python:latest'
    errorStrategy 'ignore'

    input:
    path raw

    output:
    path "${raw.baseName}-custom_headers.csv"

    """
    bruker_extract_custom_headers.py -d_folder ${raw} -out_csv ${raw.baseName}-custom_headers.csv ${params.ccff_header_in_d}
    """
}
