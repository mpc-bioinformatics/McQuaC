#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

    input:
    path raw_folder
    
    output:
    path "${raw_folder.baseName}.mzML"

    """
    tdf2mzml.py -i ${raw_folder} -o ${raw_folder.baseName}.mzML --ms1_type centroid
    """
}