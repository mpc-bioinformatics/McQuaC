/*
 * Extract XICs for spike-in/reference peptides from Bruker .d-folders, based on a
 * ThermoRawFileParser-compatible XIC extraction configuration (see XICEXTRACTIONCONFIG)
 **/
process BRUKERXICEXTRACTION {
    tag "${meta.id}"
    label 'process_medium'

    stageInMode 'copy'  // needed due to alphatims

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-7a956ec"

    input:
    tuple val(meta), path(d_folder), path(xic_config)

    output:
    tuple val(meta), path("*.json"), emit: xic
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers extract-xic-bruker \\
        ${args} \\
        -d_folder ${d_folder} \\
        -in_json ${xic_config} \\
        -out_json ${prefix}.json

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

    echo '{}' > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
        macproqc-helpers: \$(python -m macproqc_helpers --version | sed 's;__main__.py ;;')
    END_VERSIONS
    """
}
