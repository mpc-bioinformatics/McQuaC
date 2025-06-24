#!/usr/bin/env nextflow

/**
 * Workflows and processes for peptide identification using Comet
 */


nextflow.enable.dsl=2

params.identification__comet_threads = 8
// Memory per comet search
// ~6 GB for 35000 MS and 35 MB of FASTA
// Virtual and real memory were roughly equal
params.identification__comet_mem = "10 GB"
// If true, decoys are generated and searched against
params.identification__generate_decoys = true
// Method to generate decoys. OpenMS allows either: 'reverse' or 'shuffle'
params.identification__decoy_method = 'shuffle'

/*
 * Identifies peptides in MS/MS spectra using Comet
 * @param mzmls Path to mzML files
 * @param fasta_file Path to mzML file
 * @param comet_params Path to comet params file
 *
 * @return pepxmls Path to pepXML files
 */
workflow identification_with_comet {
    take:
        mzmls
        fasta_file
        mcquac_params_file
        search_labelled
        fasta_output_folder

    main:
        comet_params_file = create_fresh_comet_params()
        adj_comet_params = adjust_comet_params(comet_params_file, mcquac_params_file, search_labelled)

        if (params.identification__generate_decoys) {
            fasta_file = generate_decoy_database(fasta_file)
        }

        id_results = comet_search(mzmls, fasta_file, adj_comet_params, fasta_output_folder)
    
    emit:
        mzids = id_results.mzids
}

/**
 * Creates a fresh comet params file
 */
process create_fresh_comet_params {
    label 'comet_image'

    output:
    path "comet.params"

    script:
    """
    # create a new comet params file
    comet -p
    mv comet.params.new comet.params
    """
}

/**
 * Adjust comet.params to have the correct output files and threads are limited
 *
 * @param comet_params Comet parameter file
 * @param threads Threads to use for each comet search
 * @return adjusted comet parameter file
 */
process adjust_comet_params {
    label 'comet_image'

    input:
    path comet_params_file
    path mcquac_params_file
    val search_labelled

    output:
    path "adjusted.comet.params"

    script:
    """
    # set the number of threads
    sed -i 's/^num_threads.*/num_threads = ${params.identification__comet_threads} /' ${comet_params_file}

    adjust_comet_params.py -json_in ${mcquac_params_file} -comet_params ${comet_params_file} -params_out adjusted.comet.params -search_labelled ${search_labelled}
    """
}

/*
 * Generates a decoy database from the given FASTA file
 * @param fasta Path to fasta file
 * @param fasta_file Path to mzML file
 * @param config_file Path to comet params file
 *
 * @return pepxml Path to pepXML file
 */
process generate_decoy_database {
    label 'mcquac_image'

    input:
    path fasta

    output:
    path "${fasta.baseName}_with_decoys.fasta"

    script:
    """
    DecoyDatabase -in ${fasta} -out ${fasta.baseName}_with_decoys.fasta -method ${params.identification__decoy_method} -decoy_string DECOY_
    """
}


/*
 * Identifies peptides in MS/MS spectra using Comet
 * @param mzml Path to mzML file
 * @param fasta_file Path to mzML file
 * @param config_file Path to comet params file
 *
 * @return mzid Path to mzid file
 */
process comet_search {
    label 'comet_image'
    
    cpus { params.identification__comet_threads }
    memory { params.identification__comet_mem }

    publishDir "${fasta_output_folder}/", mode: 'copy', pattern: "*.fasta"  // Publish the FASTA file, which was used for the search

    input:
    path mzml
    path input_fasta
    path comet_params_file
    path fasta_output_folder

    output:
    path "${mzml.baseName}.mzid", emit: mzids
    path input_fasta, emit: fasta_file

    script:
    """
    comet -P${comet_params_file} -D${input_fasta} ${mzml}
    """
}
