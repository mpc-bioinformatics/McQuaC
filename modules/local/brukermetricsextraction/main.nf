/*
 * Extract specific headers from Bruker dot d folders
 **/
process BRUKERMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    stageInMode 'copy'  // needed due to pyopenms

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-7a956ec"

    input:
    tuple val(meta), path(dotd_bruker_folder)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // alphatims.bruker.TimsTOF (used in macproqc_helpers) dispatches on the folder name's
    // suffix, requiring it to end in ".d". the staged input folder may not (e.g.
    // when staged under a pipeline sample ID), so alias it to a ".d"-suffixed name.
    def orig_name = dotd_bruker_folder.name.toString()
    def dotd_name = orig_name.endsWith('.d') ? orig_name : "${orig_name}.d"
    def link_cmd = dotd_name == orig_name ? '' : "ln -s ${dotd_bruker_folder} ${dotd_name}"
    """
    ${link_cmd}
    python -m macproqc_helpers collect-metrics-from-bruker \\
        ${args} \\
        -d_folder ${dotd_name} \\
        -out_hdf5 ${prefix}.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
        macproqc-helpers: \$(python -m macproqc_helpers --version | sed 's;__main__.py ;;')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
        macproqc-helpers: \$(python -m macproqc_helpers --version | sed 's;__main__.py ;;')
    END_VERSIONS
    """
}
