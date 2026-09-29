/*
 * Collect spike-in/reference peptide metrics (maximum XIC intensity, retention time at
 * that maximum, number of identifications) from extracted XICs and identifications, and
 * store them in HDF5 format
 **/
process SPIKEINMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-7a956ec"

    input:
    tuple val(meta), path(xic_json), path(identifications)
    path spike_ins_table

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -m macproqc_helpers collect-spikein-metrics \\
        ${args} \\
        -itrfp_json ${xic_json} \\
        -iidentifications ${identifications} \\
        -ispikeins ${spike_ins_table} \\
        -ohdf5 ${prefix}.hdf5

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
