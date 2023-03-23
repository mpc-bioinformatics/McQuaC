#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.mgf_folder = "$PWD/MGFs" // Input-Directory of MGFs, which should be used for identification
params.fasta_file = "proteins.fasta" // Database (FASTA-file) used for identification with decoys (prefixed with "DECOY_") or without decoys!
params.tda = 1 // 0 --> No Target-Decoy appraoch | 1 --> Target-Decoy appraoch (comet automatically prefixes decoys with "DECOY_" )
params.search_parameter_file = "$PWD/example_configurations/comet_config.txt" //Search Parameters for Comet

// Optional Parameters
params.outdir = "$PWD/results"  // Output-Directory of the Identification Results. Here it is <Input_File>.mzid
params.num_parallel_searches = 1 //Runtime.runtime.availableProcessors()

workflow {
    // Get all MGF files which should be identified
    mgfs = Channel.fromPath(params.mgf_folder + "/*.mgf")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.search_parameter_file)

    // Combined channel replicated the indexed fasta for each MGF to be reused
    combined_channel = fasta_file
        .combine(modifications_file)
        .combine(mgfs)
        
    // Start search
    comet_search_mgf(combined_channel)
}


process comet_search_mgf {
    maxForks params.num_parallel_searches
    publishDir "${params.outdir}/idents", mode:'copy'

    input:
    tuple file(input_fasta), file(mod_file), file(mgf_file)

    output: 
    file "${mgf_file.baseName}.mzid"

    """
    sed 's/^decoy_search.*/decoy_search = ${params.tda} /' ${mod_file} > ${mod_file.baseName}_new.txt
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^decoy_prefix.*/decoy_prefix = DECOY_/' ${mod_file.baseName}_new.txt


    comet.linux_v2023.01.1.exe -P${mod_file.baseName}_new.txt -D${input_fasta} ${mgf_file}
    """  
}
