#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ctm_input_thermo_raws = "$PWD/raws"  // Folder of Thermo-RAW-files

// Optional Parameters
params.ctm_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ctm_additional_params = ""  // Addtional parameters for the TRFP
params.ctm_num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

// Standalone Workflow
workflow {
    rawfiles = Channel.fromPath(params.thermo_raws + "/*.raw")
    convert_to_mgf(rawfiles)
}


workflow convert_to_mgf_mzml {
    take:
        raw_files // a list of raw_files
    main:
        // Convert the file to MGF
        convert_raw_via_thermorawfileparser(raw_files)
    emit:
        convert_raw_via_thermorawfileparser.out[0]
        convert_raw_via_thermorawfileparser.out[1]
}

process convert_raw_via_thermorawfileparser {
    maxForks params.ctm_num_procs_conversion
    stageInMode "copy"

    publishDir "${params.ctm_outdir}/", mode:'copy'

    input:
    file raw

    output:
    file "${raw.baseName}.mgf"
    file "${raw.baseName}.mzML"

    """
    run_thermorawfileparser.sh ${params.ctm_additional_params} --format=0 --output_file=${raw.baseName}.mgf --input=${raw} 
    run_thermorawfileparser.sh ${params.ctm_additional_params} --format=1 --output_file=${raw.baseName}.mzML --input=${raw} 
    """
}
