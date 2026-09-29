/*
 * Merge all per-metric HDF5 files produced for a single run into one combined HDF5 file
 **/
process COMBINEHDF5 {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-7a956ec"

    input:
    tuple val(meta), path(hdf5_files)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers combine-hdf5 \\
        ${args} \\
        -hdf_out_name ${prefix}.hdf5 \\
        ${hdf5_files}

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
