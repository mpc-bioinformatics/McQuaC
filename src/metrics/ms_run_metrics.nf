#!/usr/bin/env nextflow
/**
 * Workflows for collecting various metadata of the MS runs from RAW files and mzML files.
 * e.g. pump preassure, ion injection time, ...
 */

nextflow.enable.dsl=2

// Set if you want to extract specific calibrants in Bruker raw measurements, otherwise the default is used 
// (622.0290 m/z, 922.009798 m/z and 1221.990637 m/z with a 10 m/z and 0.1 1/k0 tolerance).
// Have a look into the corresponding python script for the headers.


/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 * - Thermo: Lock Mass, Ion Injection Time, 
 * - Bruker: Vacuum_CurrentFore, Vacuum_CurrentHigh, etc.
 *
 * @param thermo_raw_files Channel of Thermo raw files
 * @param bruker_raw_files Channel of Bruker .d-folders
 * @return Channel of headers HDF5 files
 */
workflow get_headers {
    take:
        thermo_raw_files
        thermo_memory_limit
        thermo_headers
        bruker_raw_files
        bruker_memory_limit
        bruker_headers_to_parse

    main:
        // Extract Information directly from the RAW-file using fisher-py
        thermo_headers = extract_headers_from_thermo_raw_files(thermo_raw_files, thermo_memory_limit, thermo_headers)
        bruker_headers = extract_headers_from_bruker_raw_files(bruker_raw_files, bruker_memory_limit, bruker_headers_to_parse)
        headers = thermo_headers.concat(bruker_headers)

    emit:
        headers
}

/**
 * Get additional data from mzML, e.g. MS1_Density-quartiles, MS2_Density-quartiles, RT-TIC-quartiles, ...
 *
 * @param mzmlfiles Channel of mzML files
 * @return Channel of headers HDF5 files
 */
workflow get_mzml_infos {
    take:
        mzmlfiles
        central_params
        memory_limit

    main:
        informations = extract_data_from_mzml(mzmlfiles, central_params, memory_limit)

    emit:
        informations
}

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 *
 * @param thermo_raw_files Channel of Thermo raw files
 * @return HDF5 with extracted headers
 */
process extract_headers_from_thermo_raw_files {
    label 'mcquac_image'
    cpus 1
    memory "${memory_limit}"

    errorStrategy 'ignore'
    stageInMode 'copy'

    input:
    path raw
    val memory_limit
    val thermo_headers

    output:
    path "${raw.baseName}-custom_headers.hdf5"

    script:
    """
    # Pythonnet sometimes fails to exit and throws a mono error
    extract_thermo_headers.py -raw ${raw} -out_hdf5 ${raw.baseName}-custom_headers.hdf5 ${thermo_headers} || true

    # Fail Check if no content was written
    if ! [ -s "${raw.baseName}-custom_headers.hdf5" ];then
        rm ${raw.baseName}-custom_headers.hdf5
    fi
    """
}

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 *
 * @param bruker_raw_files Channel of Bruker .d-folders
 * @return HDF5 with extracted headers
 */
process extract_headers_from_bruker_raw_files {
    label 'mcquac_image'

    cpus 1
    memory "${memory_limit}"

    errorStrategy 'ignore'

    input:
    path raw
    val memory_limit
    val bruker_headers_to_parse

    output:
    path "${raw.baseName}-custom_headers.hdf5"

    script:
    """
    extract_bruker_headers.py -d_folder ${raw} -out_hdf5 ${raw.baseName}-custom_headers.hdf5 ${bruker_headers_to_parse} ${params.ms_run_metrics__bruker_calibrants}

    """
}

/**
 * Get additional data from mzML, e.g. MS1_Density-quartiles, MS2_Density-quartiles, RT-TIC-quartiles, ...
 *
 * @param mzml Channel of mzML files
 * @return HDF5 file with the extracted data
 */
process extract_data_from_mzml {
    label 'mcquac_image'

    cpus 1
    memory "${memory_limit}"

    input:
    path mzml
    path central_params    // JSON file with central McQuaC parameters
    val memory_limit

    output:
    path "${mzml.baseName}-mzml_info.hdf5"

    script:
    """
    extract_data_from_mzml.py -mzml ${mzml} -out_hdf5 ${mzml.baseName}-mzml_info.hdf5 -mcquac_params ${central_params}
    """
}
