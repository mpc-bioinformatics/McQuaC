
process QCVISUALIZATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'ghcr.io/mpc-bioinformatics/macproqc-helpers:latest'}"

    input:
    tuple val(meta), path(hdf5_files)
    val outputdir
    path spike_ins_table
    val figure_format
    val output_table_type
    val spikeins
    val rt_unit
    val output_column_order
    val spikein_columns
    val height_barplots
    val width_barplots
    val height_pca
    val width_pca
    val height_ionmaps
    val width_ionmaps

    output:
    tuple val(meta), path("${outputdir}/*.json"), emit: jsons, optional: true
    tuple val(meta), path("${outputdir}/*.html"), emit: htmls, optional: true
    tuple val(meta), path("${outputdir}/*.{csv,tsv,xlsx}"), emit: output_tables, optional: false
    tuple val(meta), path("${outputdir}/fig13_MS1_map"), emit: ms1_maps, optional: false
    tuple val(meta), path("${outputdir}/fig16_additional_headers"), emit: additional_plots, optional: true
    tuple val(meta), path("${outputdir}/fig17_BRUKER_calibrants"), emit: bruker_calibrants, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''   // free-form, optional extra CLI arguments
    def spike_ins_tab = spike_ins_table ? "-spike_ins_table ${spike_ins_table}" : ''
    def spikeins_arg = spikeins ? '-spikeins' : ''
    def output_column_order_arg = output_column_order ? "-output_column_order ${output_column_order}" : ''
    hdf5_files_list = hdf5_files.join(' ')
    """
    mkdir -p "${outputdir}"
    python -m macproqc_helpers visualize \\
        -hdf5_files ${hdf5_files_list} \\
        -output "${outputdir}" \\
        -figure_format ${figure_format} \\
        -output_table_type ${output_table_type} \\
        ${spikeins_arg} \\
        -RT_unit ${rt_unit} \\
        ${output_column_order_arg} \\
        -spikein_columns "${spikein_columns}" \\
        -height_barplots ${height_barplots} \\
        -width_barplots ${width_barplots} \\
        -height_pca ${height_pca} \\
        -width_pca ${width_pca} \\
        -height_ionmaps ${height_ionmaps} \\
        -width_ionmaps ${width_ionmaps} \\
        ${spike_ins_tab} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS

    """


    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch Figure1.json
    touch Outputtable.csv
    touch fig13_MS1_map
    touch fig16_additional_headers
    touch fig17_BRUKER_calibrants

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
