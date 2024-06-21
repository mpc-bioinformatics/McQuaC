#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Identifies peptides in MS/MS spectra using Comet
 * @param mzmls Path to mzML files
 * @param fasta_file Path to mzML file
 * @param config_file Path to comet params file
 *
 * @return pepxmls Path to pepXML files
 */
workflow identification_with_comet {
    take:
        mzmls
        fasta_file
        config_file

    main:
        pepxmls = comet_search(mzmls, fasta_file, config_file)
    emit:
        pepxmls
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
    container 'quay.io/medbioinf/comet-ms:v2024.01.0'
    cpus 8 // hardcode for now

    input:
    path(mzml)
    path(input_fasta)
    path(config_file)

    output: 
    path "${mzml.baseName}.pep.xml"

    """
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 0/' ${config_file.baseName}_new.txt
    sed -i 's/^output_pepxmlfile.*/output_pepxmlfile = 1/' ${config_file.baseName}_new.txt
    sed -i 's/^output_txtfile.*/output_txtfile = 0/' ${config_file.baseName}_new.txt
    sed -i 's/^num_threads.*/num_threads = 8 /' ${config_file.baseName}_new.txt

    comet -P${config_file.baseName}_new.txt -D${input_fasta} ${mzml}
    """
}
