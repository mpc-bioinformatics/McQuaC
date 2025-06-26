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
params.identification__generate_decoys = false
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
        mzmls  // Channel of mzML files, which are used for searching
        fasta_file  // A FASTA file, depending 
        comet_params  // Path to the comete parameter file
        search_labelled  // true if the search is labelled, false otherwise
        fasta_output_folder  // Folder where the FASTA file should be written. The one which is used for the search

    main:
        adj_comet_params = adjust_comet_params(comet_params, search_labelled)

        if (params.identification__generate_decoys) {
            generate_decoy_database(fasta_file)
            fasta_file = generate_decoy_database.out
        }

        id_results = comet_search(mzmls, fasta_file, adj_comet_params, fasta_output_folder)
    
    emit:
        mzids = id_results.mzids
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
    path comet_params
    val search_labelled

    output:
    path "adjusted.comet.params"

    """
    cp ${comet_params} adjusted.comet.params

    sed -i 's/^num_threads.*/num_threads = ${params.identification__comet_threads} /' adjusted.comet.params

    sed -i 's/^output_sqtfile.*/output_sqtfile = 0/' adjusted.comet.params
    sed -i 's/^output_txtfile.*/output_txtfile = 0/' adjusted.comet.params
    sed -i 's/^output_pepxmlfile.*/output_pepxmlfile = 0/' adjusted.comet.params
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 1/' adjusted.comet.params
    sed -i 's/^output_percolatorfile.*/output_percolatorfile = 0/' adjusted.comet.params

    if [ ${params.identification__decoy_method} = true ];       # hardcoded for now
    then
        sed -i 's/^decoy_search.*/decoy_search  = 0/' adjusted.comet.params
    fi
    if [ ${search_labelled} = true ];       # hardcoded for now
    then
        sed -i 's/^add_K_lysine.*/add_K_lysine  = 8.014199/' adjusted.comet.params
        sed -i 's/^add_R_arginine.*/add_R_arginine  = 10.008269/' adjusted.comet.params
    fi
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
    container { python_image }

    input:
    path fasta

    output:
    path "${fasta.baseName}_with_decoys.fasta"

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
    cpus { identification__comet_threads }
    memory { identification__comet_mem }
    publishDir "${fasta_output_folder}/", mode: 'copy', pattern: "*.fasta"  // Publish the FASTA file, which was used for the search

    input:
    path mzml
    path input_fasta
    path config_file
    path fasta_output_folder

    output:
    path "${mzml.baseName}.mzid", emit: mzids
    path input_fasta, emit: fasta_file

    """
    comet -P${config_file} -D${input_fasta} ${mzml}
    """
}
