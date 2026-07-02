process COMETCONFIG {
    tag "$meta.id"
    label 'process_single'

    publishDir path: { "${params.outdir}/comet" }, mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.14'
        : 'biocontainers/python:3.14'}"

    input:
    tuple val(meta), path(comet_config_template)

    output:
    tuple val(meta), path("${params_out_file}"), emit: params
    path "versions.yml"                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    params_out_file = "${prefix}.comet.params"
    template 'adjust_comet_params.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    params_out_file = "${prefix}.comet.params"
    """
    touch ${params_out_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
