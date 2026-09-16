/*
 * Collect spike-in/reference peptide metrics (maximum XIC intensity, retention time at
 * that maximum, number of identifications) from extracted XICs and identifications, and
 * store them in HDF5 format
 **/
process SPIKEINMETRICSEXTRACTION {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-914105c"

    input:
    tuple val(meta), path(xic_json), path(identifications)
    path spike_ins_table

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    tuple val("${task.process}"), val('macproqc-helpers'), val("sha-914105c"), topic: versions, emit: versions_spikeinmetricsextraction

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
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.hdf5
    """
}
