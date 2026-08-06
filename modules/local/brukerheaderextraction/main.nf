/*
 * Extract specific headers from Bruker dot d folders
 **/
process BRUKERHEADEREXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    stageInMode 'copy'  // needed due to pyopenms

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-914105c"

    input:
    tuple val(meta), path(dotd_bruker_folder)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc-helpers'), val("sha-914105c"), emit: versions_brukerheaderextraction

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers collect-metrics-from-bruker \\
        ${args} \\
        -d_folder ${dotd_bruker_folder} \\
        -out_hdf5 ${prefix}.hdf5
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}.hdf5
    """
}
