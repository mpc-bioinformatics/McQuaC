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
    container 'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses'
    // by mounting the parent directory we can use a symlink to the raw file in the workdir
    containerOptions { "-v ${raw_file.getParent()}:/data" }

    input:
    path raw_file

    output:
    path "${raw_file.baseName}.mzML"

    """
    wine msconvert ${raw_file} --mzML --zlib --filter "peakPicking true 1-"
    """
}


process convert_bruker_raw_files {
    container 'mfreitas/tdf2mzml'
    containerOptions { "-v ${raw_file.getParent()}:/data" }

    maxForks params.ctm_num_procs_conversion

    publishDir "${params.ctm_outdir}/", mode:'copy'

    input:
    file raw_file
    
    output:
    file("${raw_file.baseName}.mzML")

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