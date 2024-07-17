#!/usr/bin/env nextflow
nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'

/*
 * Writes out processing information
 * 
 * @param raw_files the raw files of the complete project
 * @param processed_basenames the basenames of the successfully processed files
 *
 * @return collected_information a file containing human readable information about whether teh files were processed successfully
 */
workflow output_processing_success {
    take:
    raw_files
    processed_basenames

    main:
    raw_map = raw_files.map{file -> tuple(file.baseName, file)}
	processed_ok = processed_basenames.map{name -> tuple(name, "successful")}

    info_files = write_processing_information(raw_map.join(processed_ok, remainder: true))

    collected_information = info_files.collectFile(name: "processing_info.txt", storeDir: "${params.main_outdir}")

    emit:
    collected_information
}


/**
 * Writes for each raw file whether it was processed successfully or not.  
 *
 * @param tuple of val(file_base), path(raw_file), val(information), while information is "successful" or null
 *
 * @return a file containing human readable info about the processing state
 */
process write_processing_information {
    container { python_image }

	input:
	tuple val(file_base), path(raw_file), val(information)

    output:
    path "${file_base}_processing_info.txt"

    """
    if [ "${information}" = "null" ]; then
        msg="ERROR processing"
    else
        msg="successfully processed"
    fi
    
    echo "\${msg} ${file_base} (${raw_file})" > ${file_base}_processing_info.txt
    """
}