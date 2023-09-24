#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.mgf_folder = "$PWD/MGFs" // Input-Directory of MGFs, which should be used for identification
params.fasta_file = "proteins.fasta" // Database (FASTA-file) used for identification with decoys (prefixed with "DECOY_") or without decoys!
params.tda = 1 // 0 --> No Target-Decoy appraoch | 1 --> Target-Decoy appraoch (we hardcoded the prefix "DECOY_" )
params.search_parameter_file = "$PWD/example_configurations/msgfplus_config.txt" //Search Parameters for MSGFPlus

// Optional Parameters
params.outdir = "$PWD/results"  // Output-Directory of the Identification Results. Here it is <Input_File>.mzid
params.num_parallel_searches = Runtime.runtime.availableProcessors()

// Parameters for the JVM
params.jvm_params = "Xmx3500M"

workflow {
    // Get all mgf files which should be identified
    mgfs = Channel.fromPath(params.mgf_folder + "/*.mgf")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.search_parameter_file)

    // Build indexed fasta for MSGFPLUS
    msgfplus_buildsa(fasta_file)

    // Combined channel replicated the indexed fasta for each MGF to be reused
    combined_channel = fasta_file
        .combine(modifications_file)
        .combine(mgfs)
        .combine(msgfplus_buildsa.out.toList())
        
    // Start search
    msgfplus_search_mgf(combined_channel)
}

process msgfplus_buildsa {

    input:
    file input_fasta

    output:
    file "${input_fasta.baseName}*"
    // file "${input_fasta.baseName}*.fasta"
    // file "${input_fasta.baseName}*.canno"
    // file "${input_fasta.baseName}*.cnlcp"
    // file "${input_fasta.baseName}*.csarr"
    // file "${input_fasta.baseName}*.cseq"    

    """
    buildsa_msgfplus.sh ${params.jvm_params} -d ${input_fasta} -tda ${params.tda} -decoy DECOY_
    """
}

process msgfplus_search_mgf {
    maxForks params.num_parallel_searches
    stageInMode "copy"
    publishDir "${params.outdir}/mzid", mode:'copy'

    input:
    tuple file(input_fasta), file(mod_file), file(mgf_file), file(remainder)

    output: 
    file "${mgf_file.baseName}.mzid"

    """
    touch -m ${input_fasta.baseName}*

    java -${params.jvm_params} -jar \$(get_cur_bin_dir.sh)/MSGFPLUS_v20220418/MSGFPlus.jar -d ${input_fasta} -s ${mgf_file} -tda ${params.tda} -decoy DECOY_ -conf ${mod_file} -o ${mgf_file.baseName}.mzid

    """  
}

//    # Since we copy the input files, and MS-GF+ does not like symlinks, we later remove them manually to save memory
//    rm ${input_fasta.baseName}*
//    rm ${mgf_file}
//    rm ${mod_file}