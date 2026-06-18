process THERMOHEADEREXTRACTION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"  // TODO DL it is not needed, can we remove this?
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/mpc-bioinformatics/mcquac:0.1.0':
        'ghcr.io/mpc-bioinformatics/mcquac:0.1.0' }"

    input:
    tuple val(meta), path(raw_thermo_file)
    val(thermo_extraction_headers)

    // TODO DL how to do versioning?
    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    path("versions.yml"), emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile_basename = prefix
    template 'thermo_header_extraction.py'

    stub:
    // TODO DL how to do versioning? Does it need to be in stub?
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    echo ${raw_thermo_file}
    echo ${thermo_extraction_headers}
    
    touch ${prefix}.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
        thermo_header_extraction: "0.0.1"
    END_VERSIONS
    """
}

