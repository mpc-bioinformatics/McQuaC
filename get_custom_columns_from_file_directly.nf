#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ccff_input_thermo = "$PWD/raws"  // Folder of Thermo-RAW-files

// Optional Parameters
params.ccff_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ccff_header_in_raws = '-htp "Ion Injection Time" -htp "Number of Lock Masses" -htp "Lock Mass #1" -htp "Lock Mass #2" -htp "Lock Mass #3" -htp "LM Search Window" -htp "LM Search Window" -htp "Number of LM Found" -htp "Last Locking" -htp "LM m/z-Correction"'  // Headers which we want to extract from thermo

params.ccff_header_in_raws_names = '-cn "Ion_Injection_Time" -cn "Number_of_Lock_Masses" -cn "Lock_Mass_1" -cn "Lock_Mass_2" -cn "Lock_Mass_3" -cn "LM_Search_Window" -cn "LM_Search_Window" -cn "Number_of_LM_Found" -cn "Last_Locking" -cn "LM_m_z_Correction"'  // Headers column names which we set for the output csv

// Standalone Workflow
workflow {
    rawfiles = Channel.fromPath(params.ccff_input_thermo + "/*.raw")
    get_custom_headers(rawfiles)
}


workflow get_custom_headers {
    take:
        raw_files // a list of raw_files
    main:
        // Convert the file to MGF
        retrieve_custom_headers_thermo(raw_files)
    emit:
        retrieve_custom_headers_thermo.out
}

process retrieve_custom_headers_thermo {
    maxForks 2
    publishDir "${params.ccff_outdir}/", mode:'copy'
    debug true

    input:
    file raw

    output:
    file "${raw.baseName}_____customs.csv"

    """
    thermo_extract_cutsom_headers.py -raw ${raw} ${params.ccff_header_in_raws} ${params.ccff_header_in_raws_names} -out_csv ${raw.baseName}_____customs.csv 
    """
}
