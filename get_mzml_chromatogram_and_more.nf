#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.gmc_thermo_raws = "$PWD/raws"  // Folder of Thermo-RAW-files

// Optional Parameters
params.gmc_outdir = "$PWD/results"  // Output-Directory of the mzMLs. Here it is <Input_file>.mzML
params.gmc_num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

// Standalone Workflow
workflow {
    rawfiles = Channel.fromPath(params.gmc_thermo_raws + "/*.raw")
    get_various_mzml_infos(raw_files)
}


workflow get_various_mzml_infos {
    take:
        raw_files
    main:
        // Convert the file to mzML, where MS1 (peak picked) (for feature finding)
        convert_raw_via_thermorawfileparser(raw_files)

        // Retrieve the information from mzml (using pyopenms in a python script)
        retrieve_data_from_mzml(convert_raw_via_thermorawfileparser.out)
    emit: 
        convert_raw_via_thermorawfileparser.out
        retrieve_data_from_mzml.out
}

process convert_raw_via_thermorawfileparser {
    maxForks params.gmc_num_procs_conversion
    stageInMode "copy"

    publishDir "${params.gmc_outdir}/", mode:'copy'

    input:
    file(raw)

    output:
    file("${raw.baseName}.mzML")

    """
    run_thermorawfileparser.sh --format=1 --output_file=${raw.baseName}.mzML --input=${raw} 
    rm ${raw}
    """
}

process retrieve_data_from_mzml {
    stageInMode "copy"

    publishDir "${params.gmc_outdir}/", mode:'copy'

    input:
    file(mzml)

    output:
    file("${mzml.baseName}_____mzml_info.csv")

    """
    extract_data_from_mzml.py -mzml ${mzml} -out_csv ${mzml.baseName}_____mzml_info.csv
    """
}
