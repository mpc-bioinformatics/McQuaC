/*
 * Extract specific metrics from mzML files
 **/
process MZMLMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_low'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-914105c"

    input:
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc-helpers'), val("sha-914105c"), emit: versions_mzmlmetricsextraction

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