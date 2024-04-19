#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ic_mzmls = "$PWD/MZMLs" // Input-Directory of mzmls, which should be used for identification
params.ic_fasta_file = "proteins.fasta" // Database (FASTA-file) used for identification with decoys (prefixed with "DECOY_") or without decoys!
params.ic_search_parameter_file = "${baseDir}/example_configurations/high-high.comet.params" //Search Parameters for Comet

// Optional Parameters
params.ic_tda = 1 // 0 --> No Target-Decoy appraoch | 1 --> Target-Decoy appraoch (comet automatically prefixes decoys with "DECOY_" )
params.ic_outdir = "$PWD/results"  // Output-Directory of the Identification Results. Here it is <Input_File>.mzid
params.ic_num_parallel_threads_per_search = 4


workflow {
    // Get all MGF files which should be identified
    mzmls = Channel.fromPath(params.ic_mzmls + "/*.mzML")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.ic_fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.ic_search_parameter_file)

    ident_via_comet(mzmls, fasta_file, modifications_file)
}

workflow ident_via_comet {
    take:
        mzmls
        fasta_file
        config_file

    main:
        // Start search
        idxmls = comet_search(mzmls, fasta_file, config_file)
    emit:
        idxmls
}

process comet_search {
    container 'mpc/nextqcflow-comet:latest'

    cpus params.ic_num_parallel_threads_per_search
    publishDir "${params.ic_outdir}/idents", mode:'copy'

    input:
    path(mzml)
    path(input_fasta)
    path(config_file)


    output: 
    path "${mzml.baseName}.pep.xml"

    """
    sed 's/^decoy_search.*/decoy_search = ${params.ic_tda} /' ${config_file} > ${config_file.baseName}_new.txt
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 0/' ${config_file.baseName}_new.txt
    sed -i 's/^decoy_prefix.*/decoy_prefix = DECOY_/' ${config_file.baseName}_new.txt
    sed -i 's/^database_name.*/database_name = ${input_fasta}/' ${config_file.baseName}_new.txt
    sed -i 's/^output_pepxmlfile.*/output_pepxmlfile = 1/' ${config_file.baseName}_new.txt
    sed -i 's/^output_txtfile.*/output_txtfile = 0/' ${config_file.baseName}_new.txt
    sed -i 's/^num_threads.*/num_threads = ${params.ic_num_parallel_threads_per_search} /' ${config_file.baseName}_new.txt


    # We run the Comet Adapter since PIA does not work with comet output and need to mimic everything to satisfy OpenMS
    # DecoyDatabase -in ${input_fasta} -out ${input_fasta.baseName}_rev.fasta -decoy_string "DECOY_" -method "reverse"
    # CometAdapter -PeptideIndexing:unmatched_action warn  -in ${mzml} -out ${mzml.baseName}.idXML -database ${input_fasta.baseName}_rev.fasta -comet_executable ${workflow.projectDir}/bin/comet.linux_v2022.01.0.exe -default_params_file ${config_file.baseName}_new.txt

    comet -P${config_file.baseName}_new.txt -D${input_fasta} ${mzml}
    """
}
