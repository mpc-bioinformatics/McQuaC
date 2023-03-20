#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.thermo_raws = "$PWD/raws"  // Folder of Thermo-RAW-files

// Optional Parameters

// Output Directory
params.outdir = "$PWD/results"  // Output-Directory of the mzMLs. Here it is <Input_file>.mzML

params.num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

workflow {
    // Convert the file to mzML, where MS1 (peak picked) (for feature finding)
    rawfiles = Channel.fromPath(params.thermo_raws + "/*.raw")
    convert_raw_via_thermorawfileparser(rawfiles)

    // Retrieve the information from mzml (using pyopenms in a python script)
    retrieve_data_from_mzml(convert_raw_via_thermorawfileparser.out)
}

process convert_raw_via_thermorawfileparser {
    maxForks params.num_procs_conversion
    stageInMode "copy"

    publishDir "${params.outdir}/", mode:'copy'

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
    maxForks 1  // TODO Please remove limit, my laptop is not that powerful! 
    stageInMode "copy"

    publishDir "${params.outdir}/", mode:'copy'

    input:
    file(mzml)

    output:
    file("${mzml.baseName}_mzml_info.csv")

    """
    extract_data_from_mzml.py -mzml ${mzml} -out_csv ${mzml.baseName}_mzml_info.csv
    """
}
