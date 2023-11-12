#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ccff_input_thermo = "$PWD/raws"  // Folder of Thermo-RAW-files

// Optional Parameters
params.ccff_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ccff_header_in_raws = ""  // Headers which we want to extract from thermo

params.ccff_header_in_d = "-htp 'Vacuum_CurrentFore' -htp 'Vacuum_Extra4thGauge' -htp 'Vacuum_CurrentHigh' -htp 'Vacuum_CurrentFunnel' -htp 'Digitizer_CurrentTemp' -htp 'TOF_DeviceTempCurrentValue1' -htp 'TOF_DeviceTempCurrentValue2'"
params.ccff_header_in_d_names = "-cn 'BRUKER_Vacuum_CurrentFore_pickle_zlib' -cn 'BRUKER_Vacuum_Extra4thGauge_pickle_zlib' -cn 'BRUKER_Vacuum_CurrentHigh_pickle_zlib' -cn 'BRUKER_Vacuum_CurrentFunnel_pickle_zlib' -cn 'BRUKER_Digitizer_CurrentTemp_pickle_zlib' -cn 'BRUKER_TOF_DeviceTempCurrentValue1_pickle_zlib' -cn 'BRUKER_TOF_DeviceTempCurrentValue2_pickle_zlib'"

// Standalone Workflow
workflow {
    rawfiles = Channel.fromPath(params.ccff_input_thermo + "/*.raw")
    get_custom_headers(rawfiles)
}


workflow get_custom_headers {
    take:
        raw_files // a list of raw_files
    main:
        // Extract Information directly from the RAW-file using fisher-py
        retrieve_custom_headers(raw_files)
    emit:
        retrieve_custom_headers.out
}

process retrieve_custom_headers {
    stageInMode "copy"
    publishDir "${params.ccff_outdir}/", mode:'copy'

    input:
    file raw

    output:
    file "${raw.baseName}_____customs.csv"

    """
    if [[ "${raw}" == *.raw ]]; then
        # Pythonnet sometimes fails to exit and throws a mono error
        thermo_extract_custom_headers.py -raw ${raw} ${params.ccff_header_in_raws} -out_csv ${raw.baseName}_____customs.csv || true

        # Fail Check if no content was written
        if ! [ -s "${raw.baseName}_____customs.csv" ];then
            rm ${raw.baseName}_____customs.csv
        fi
    fi

    if [[ "${raw}" == *.d ]]; then
        bruker_extract_custom_headers.py -d_folder ${raw} -out_csv ${raw.baseName}_____customs.csv ${params.ccff_header_in_d} ${params.ccff_header_in_d_names}
    fi 
 
    """
}
