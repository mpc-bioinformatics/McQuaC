#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ctm_input_spectra = "$PWD/raws"  // Folder of Thermo (.RAW) or Bruker (.d) Spectra files/folders

// Optional Parameters
params.ctm_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ctm_additional_params = ""  // Addtional parameters for the TRFP
params.ctm_num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

// Standalone Workflow
workflow {
    rawfiles = Channel.fromPath(params.ctm_input_spectra + "/*.{raw,d}", type: "any")
    convert_to_mgf_mzml(rawfiles)
}


workflow convert_to_mgf_mzml {
    take:
        raw_files // a list of raw_files
    main:
        // Convert the file to MGF
        convert_spectra(raw_files)
    emit:
        convert_spectra.out[0]
        convert_spectra.out[1]
}

process convert_spectra {
    maxForks params.ctm_num_procs_conversion
    stageInMode "copy"

    publishDir "${params.ctm_outdir}/", mode:'copy'

    input:
    file raw_spectra

    output:
    file("${raw_spectra.baseName}.mgf")
    file("${raw_spectra.baseName}.mzML")

    """
    if [[ "${raw_spectra}" == *.raw ]]; then
        thermorawfileparser ${params.ctm_additional_params} --format=0 --output_file=${raw_spectra.baseName}.mgf --input=${raw_spectra} 
        thermorawfileparser ${params.ctm_additional_params} --format=1 --output_file=${raw_spectra.baseName}.mzML --input=${raw_spectra} 
    fi

    if [[ "${raw_spectra}" == *.d ]]; then
        tdf2mzml_wrapper -i ${raw_spectra} -o ${raw_spectra.baseName}.mzML --ms1_type centroid
        
        alphatims export mgf -o . ${raw_spectra}

        # TODO We need to find a better way. MGF is optional. This is not needed for a DIA analysis
        touch ${raw_spectra.baseName}.mgf
    fi 
    """
}

