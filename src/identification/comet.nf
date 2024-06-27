#!/usr/bin/env nextflow
nextflow.enable.dsl=2

comet_image = 'quay.io/medbioinf/comet-ms:v2024.01.0'
comet_threads = 8 // hardcoded for now

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
        comet_params
        search_labelled

    main:
        adj_comet_params = adjust_comet_params(comet_params, search_labelled)
        pepxmls = comet_search(mzmls, fasta_file, adj_comet_params)
    
    emit:
        pepxmls
}


/**
 * Adjust comet.params to have the correct output files and threads are limited
 *
 * @param comet_params Comet parameter file
 * @param threads Threads to use for each comet search
 * @return adjusted comet parameter file
 */
process adjust_comet_params {
    container { comet_image }

    input:
    path comet_params
    val search_labelled

    output:
    path "adjusted.comet.params"

    """
    cp ${comet_params} adjusted.comet.params

    sed -i 's/^num_threads.*/num_threads = ${comet_threads} /' adjusted.comet.params

    sed -i 's/^output_sqtfile.*/output_sqtfile = 0/' adjusted.comet.params
    sed -i 's/^output_txtfile.*/output_txtfile = 0/' adjusted.comet.params
    sed -i 's/^output_pepxmlfile.*/output_pepxmlfile = 1/' adjusted.comet.params
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 0/' adjusted.comet.params
    sed -i 's/^output_percolatorfile.*/output_percolatorfile = 0/' adjusted.comet.params

    if [ ${search_labelled} = true ];       # hordcoded for now
    then
        sed -i 's/^add_K_lysine.*/add_K_lysine  = 8.014199/' adjusted.comet.params
        sed -i 's/^add_R_arginine.*/add_R_arginine  = 10.008269/' adjusted.comet.params
    fi
    """
}

/*
 * Identifies peptides in MS/MS spectra using Comet
 * @param mzml Path to mzML file
 * @param fasta_file Path to mzML file
 * @param config_file Path to comet params file
 *
 * @return pepxml Path to pepXML file
 */
process comet_search {
    container { comet_image }
    cpus { comet_threads } // hardcode for now

    input:
    path mzml
    path input_fasta
    path config_file

    output: 
    path "${mzml.baseName}.pep.xml"

    """
    comet -P${config_file} -D${input_fasta} ${mzml}
    """
}
