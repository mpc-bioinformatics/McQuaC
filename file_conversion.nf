#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ctm_input_spectra = "$PWD/raws"  // Folder of Thermo (.RAW) or Bruker (.d) Spectra files/folders

// Optional Parameters
params.ctm_outdir = "$PWD/results"  // Output-Directory of the MGFs. Here it is <Input_file>.mgf
params.ctm_additional_params = ""  // Addtional parameters for the TRFP
// Memory for the Thermo Raw File Parser, used 24 GB for a Raw file with 257409 MS scans 
// and 4GB for a Raw file with 11352 MS scans (measured with `/usr/bin/time -v ...`). 10 GB seems legit for most cases.
// Based on max virtual memory 
params.file_conversion__thermo_raw_conversion_mem = "10 GB"
// Memory for the tdf2mzml, used 0.39 GB for a Raw file with 298748 MS scans 
/// and 0.14GB for a Raw file with 35023 MS scans (measured with `/usr/bin/time -v ...`). 5 GB seems more then enough.
// Based on max virtual memory 
params.file_conversion__bruker_raw_conversion_mem = "5 GB"

// Standalone Workflow
workflow {
    thermo_raw_files = Channel.fromPath(params.ctm_input_spectra + "/*.raw")
    bruker_raw_files = Channel.fromPath(params.ctm_input_spectra + "/*.d", type: "dir")
    convert_mzml(thermo_raw_files, bruker_raw_files)
}

/**
 * Convert raw files (Thermo Fisher .raw-files and Bruker tdf-files) to mzML files
 */
workflow convert_to_mzml {
    take:
        thermo_raw_files // a list of raw_files
        bruker_raw_files // a list of raw_files
    main:
        thermo_mzmls = convert_thermo_raw_files(thermo_raw_files)
        bruker_mzmls = convert_bruker_raw_files(bruker_raw_files)
        mzmls = thermo_mzmls.concat(bruker_mzmls)
    emit:
        mzmls
}

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


process convert_thermo_raw_files {
    container 'quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0'
    // Thermo Raw File parser is currently limited to 2 CPUs, see:
    // * https://github.com/compomics/ThermoRawFileParser/issues/23
    // * https://github.com/compomics/ThermoRawFileParser/issues/95
    cpus 2
    memory params.file_conversion__thermo_raw_conversion_mem

    input:
    path raw_file

    output:
    path "${raw_file.baseName}.mzML"

    """
    thermorawfileparser ${params.ctm_additional_params} --format=2 --output_file=${raw_file.baseName}.mzML --input=${raw_file}
    """
}


process convert_bruker_raw_files {
    container 'mfreitas/tdf2mzml'
    containerOptions { "-v ${raw_file.getParent()}:/data" }
    // Uses all cores
    cpus Runtime.runtime.availableProcessors()
    memory params.file_conversion__bruker_raw_conversion_mem


    publishDir "${params.ctm_outdir}/", mode:'copy'

    input:
    path raw_file
    
    output:
    path "${raw_file.baseName}.mzML"

    """
    tdf2mzml.py -i ${raw_file} -o ${raw_file.baseName}.mzML --ms1_type centroid
    """
}

/**
 * Convert mzML files to MGF files
 * TODO: Is peak picking necessary? Only when mzmls not peak picked?
 */
process convert_mzml {
    container 'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses'
    // by mounting the parent directory we can use a symlink to the raw file in the workdir
    containerOptions { "-v ${raw_file.getParent()}:/data" }

    input:
    path mzml

    output:
    path "${raw_spectra.baseName}.mgf"

    """
    wine msconvert ${mzml} --mgf --filter "peakPicking true 1-"
    """
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