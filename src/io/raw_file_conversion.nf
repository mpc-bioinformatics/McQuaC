#!/usr/bin/env nextflow

/**
 * Workflows and processes file conversions, e.g. from vendor file formats to open standards like mzML 
 */

nextflow.enable.dsl=2

// Memory for the Thermo Raw File Parser, used 24 GB for a Raw file with 257409 MS scans 
// and 4GB for a Raw file with 11352 MS scans (measured with `/usr/bin/time -v ...`). 10 GB seems legit for most cases.
// Based on max virtual memory 
params.file_conversion__thermo_raw_conversion_mem = "10 GB"
// Memory for the tdf2mzml, used 0.39 GB for a Raw file with 298748 MS scans 
// and 0.14GB for a Raw file with 35023 MS scans (measured with `/usr/bin/time -v ...`). 5 GB seems more then enough.
// Based on max virtual memory 
params.file_conversion__bruker_raw_conversion_mem = "5 GB"

/**
 * Convert raw files (Thermo Fisher .raw-files and Bruker tdf-files) to mzML files
 * @params thermo_raw_files Thermo Fisher .raw-files
 * @params bruker_raw_folders Bruker .d-folders
 */
workflow convert_raws_to_mzml {
    take:
        thermo_raw_files
        bruker_raw_folders
    main:
        thermo_mzmls = convert_thermo_raw_files(thermo_raw_files)
        bruker_mzmls = convert_bruker_raw_folders(bruker_raw_folders)
        mzmls = thermo_mzmls.concat(bruker_mzmls)
    emit:
        mzmls
}

/**
 * Convert raw file (Thermo Fisher .raw-files) to mzML files
 * @params thermo_raw_files A list of Thermo Fisher .raw-files
 *
 * @return mzML files
 */
process convert_thermo_raw_files {
    container 'quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0'
    errorStrategy 'ignore'
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
    thermorawfileparser --format=2 --output_file=${raw_file.baseName}.mzML --input=${raw_file}
    """
}

/**
 * Convert raw files (Bruker .d-folder) to mzML files
 * @params raw_folder Bruker .d-folder
 *
 * @return mzML file
 */
process convert_bruker_raw_folders {
    container 'mfreitas/tdf2mzml'
    containerOptions { "-v ${raw_folder.getParent()}:/data" }
    errorStrategy 'ignore'
    // Uses all cores
    cpus 4
    memory params.file_conversion__bruker_raw_conversion_mem

    input:
    path raw_folder
    
    output:
    path "${raw_folder.baseName}.mzML"

    """
    tdf2mzml.py -i ${raw_folder} -o ${raw_folder.baseName}.mzML --ms1_type centroid
    """
}
