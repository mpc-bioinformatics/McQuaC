#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.ic_raw_folder = "$PWD/MZMLs" // Input-Directory of MGFs, which should be used for identification
params.ic_fasta_file = "proteins.fasta" // Database (FASTA-file) used for identification with decoys (prefixed with "DECOY_") or without decoys!
params.ic_search_parameter_file = "$PWD/example_configurations/comet_config.txt" //Search Parameters for Comet

// Optional Parameters
params.ic_tda = 1 // 0 --> No Target-Decoy appraoch | 1 --> Target-Decoy appraoch (comet automatically prefixes decoys with "DECOY_" )
params.ic_outdir = "$PWD/results"  // Output-Directory of the Identification Results. Here it is <Input_File>.mzid
params.ic_num_parallel_threads_per_search = 4


workflow {
    // Get all MGF files which should be identified
    mgfs = Channel.fromPath(params.ic_raw_folder + "/*.mgf")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.ic_fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.ic_search_parameter_file)

    ident_via_comet(mgfs, fasta_file, modifications_file)
}

workflow ident_via_comet {
    take:
        mgfs
        fasta_file
        modification_file

    main:
        // Combined channel replicated the indexed fasta for each MGF to be reused
        combined_channel = fasta_file
            .combine(modification_file)
            .combine(mgfs)
        
        // Start search
        comet_search(combined_channel)
    emit:
        comet_search.out
}

process comet_search {
    cpus params.ic_num_parallel_threads_per_search
    publishDir "${params.ic_outdir}/idents", mode:'copy'
    stageInMode "copy"

    input:
    tuple file(input_fasta), file(mod_file), file(mzml)

    output: 
    file "${mzml.baseName}.idXML"

    """
    sed 's/^decoy_search.*/decoy_search = ${params.ic_tda} /' ${mod_file} > ${mod_file.baseName}_new.txt
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^decoy_prefix.*/decoy_prefix = DECOY_/' ${mod_file.baseName}_new.txt
    sed -i 's/^database_name.*/database_name = ${input_fasta}/' ${mod_file.baseName}_new.txt
    sed -i 's/^output_pepxmlfile.*/output_pepxmlfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^output_txtfile.*/output_txtfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^num_threads.*/num_threads = ${params.ic_num_parallel_threads_per_search} /' ${mod_file.baseName}_new.txt


    # We run the Comet Adapter since PIA does not work with comet output and need to mimic everything to satisfy OpenMS
    #run_decoydatabase.sh -in ${input_fasta} -out ${input_fasta.baseName}_rev.fasta -decoy_string "DECOY_" -method "reverse"
    #run_cometadapter.sh -PeptideIndexing:unmatched_action warn  -in ${mzml} -out ${mzml.baseName}.idXML -database ${input_fasta.baseName}_rev.fasta -comet_executable ${workflow.projectDir}/bin/comet.linux_v2022.01.0.exe -default_params_file ${mod_file.baseName}_new.txt

    comet.linux_v2022.01.2.exe -P${mod_file.baseName}_new.txt -D${input_fasta} ${mzml}
    run_idfileconverter.sh -in ${mzml.baseName}.pep.xml -out ${mzml.baseName}.idXML 
    """  
}
