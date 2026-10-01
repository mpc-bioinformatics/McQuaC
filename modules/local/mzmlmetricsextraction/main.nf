/*
 * Extract specific metrics from mzML files
 **/
process MZMLMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-604eac1"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc_helpers'), eval('macproqc-helpers --version | sed "s;^macproqc-helpers ;;"'), topic: versions, emit: versions_macproqc_helpers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python -m macproqc_helpers collect-metrics-from-mzml \\
        ${args} \\
        -mzml ${mzml_file} \\
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
