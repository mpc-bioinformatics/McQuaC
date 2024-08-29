#!/usr/bin/env nextflow
/**
 * Workflows for collecting various metadata of the MS runs from RAW files and mzML files.
 * e.g. pump preassure, ion injection time, ...
 */

nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'

// Set if you want to extract specific headers from Thermo raw measurements, otherwise the default is used.
// Have a look into the corresponding python script for the headers.
params.ms_run_metrics__thermo_headers = ""
 // Set if you want to extract specific headers from Bruker raw measurements, otherwise the default is used.
 // Have a look into the corresponding python script for the headers.
params.ms_run_metrics__bruker_headers = ""
// Memory for the Thermo Raw File Parser, used 24 GB for a Raw file with 257409 MS scans 
/// and 4GB for a Raw file with 11352 MS scans (measured with `/usr/bin/time -v ...`). 10 GB seems legit for most cases.
params.ms_run_metrics__thermo_raw_mem = "10 GB"
// Memory for the tdf2mzml, used 0.39 GB for a Raw file with 298748 MS scans 
/// and 0.14GB for a Raw file with 35023 MS scans (measured with `/usr/bin/time -v ...`). 5 GB seems more then enough.
params.ms_run_metrics__bruker_raw_mem = "1 GB"
/// Tracing showed up to 4.7 GB virtual memory for 30000 MS scans
params.ms_run_metrics__mzml_mem = "10 GB"

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 * - Thermo: Lock Mass, Ion Injection Time, 
 * - Bruker: Vacuum_CurrentFore, Vacuum_CurrentHigh, etc.
 *
 * @param thermo_raw_files Channel of Thermo raw files
 * @param bruker_raw_files Channel of Bruker .d-folders
 * @return Channel of headers CSV files
 */
workflow get_headers {
    take:
        thermo_raw_files
        bruker_raw_files
    main:
        // Extract Information directly from the RAW-file using fisher-py
        thermo_headers = extract_headers_from_thermo_raw_files(thermo_raw_files)
        bruker_headers = extract_headers_from_bruker_raw_files(bruker_raw_files)
        headers = thermo_headers.concat(bruker_headers)
    emit:
        headers
}

/**
 * Get additional data from mzML, e.g. MS1_Density-quartiles, MS2_Density-quartiles, RT-TIC-quartiles, ...
 *
 * @param mzmlfiles Channel of mzML files
 * @return Channel of headers CSV files
 */
workflow get_mzml_infos {
    take:
        mzmlfiles
    main:
        informations = extract_data_from_mzml(mzmlfiles)
    emit:
        informations
}

/**
 * Get metadata headers from Thermo and Bruker raw files, like 
 *
 * @param thermo_raw_files Channel of Thermo raw files
 * @return CSV with 
 */
process extract_headers_from_thermo_raw_files {
    container { python_image}
    errorStrategy 'ignore'

    memory params.ms_run_metrics__thermo_raw_mem

    stageInMode 'copy'

    input:
    path raw

    output:
    path "${raw.baseName}-custom_headers.hdf5"

    """
    # Pythonnet sometimes fails to exit and throws a mono error
    extract_thermo_headers.py -raw ${raw} ${params.ms_run_metrics__thermo_headers} -out_hdf5 ${raw.baseName}-custom_headers.hdf5 || true

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
 * @return CSV with 
 */
process extract_headers_from_bruker_raw_files {
    container { python_image}
    errorStrategy 'ignore'

    memory params.ms_run_metrics__bruker_raw_mem

    input:
    path raw

    output:
    path "${raw.baseName}-custom_headers.csv"

    """
    extract_bruker_headers.py -d_folder ${raw} -out_csv ${raw.baseName}-custom_headers.csv ${params.ms_run_metrics__bruker_headers}
    """
}

/**
 * Get additional data from mzML, e.g. MS1_Density-quartiles, MS2_Density-quartiles, RT-TIC-quartiles, ...
 *
 * @param mzml Channel of mzML files
 * @return CSV file with the extracted data
 */
process extract_data_from_mzml {
    container { python_image}

    cpus 1
    memory params.ms_run_metrics__mzml_mem

    input:
    path(mzml)

    output:
    path("${mzml.baseName}-mzml_info.hdf5")

    """
    extract_data_from_mzml.py -mzml ${mzml} -out_hdf5 ${mzml.baseName}-mzml_info.hdf5
    """
}
