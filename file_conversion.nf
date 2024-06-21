#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ctm_input_spectra = "$PWD/raws"  // Folder of Thermo (.RAW) or Bruker (.d) Spectra files/folders

// Optional Parameters
params.ctm_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ctm_additional_params = ""  // Addtional parameters for the TRFP
params.ctm_num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

/**
 * Convert mzML files to MGF files
 */
workflow convert_to_mgf{
    take:
        mzmls
    main:
        mgfs = convert_mzml(mzmls)
    emit:
        mgfs
}

/**
 * Convert pepxml files to idXML files
 */
workflow convert_to_idxml{
    take:
        pepxml
    main:
        idxml = convert_pepxml_to_idxml(pepxml)
    emit:
        idxml
}

process convert_pepxml_to_idxml {
    container 'mpc/nextqcflow-python:latest'

    input:
    path pepxml

    output:
    path "${pepxml.getSimpleName()}.idXML"

    """
    IDFileConverter -in ${pepxml} -out ${pepxml.getSimpleName()}.idXML 
    """
}