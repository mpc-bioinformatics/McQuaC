
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
    def args = task.ext.args ?: ''
    def spikein_columns = task.ext.spikein_columns  // because it contains whitespace
    def spike_ins_tab = spike_ins_table ? "-spike_ins_table ${spike_ins_table}" : ''
    hdf5_files_list = hdf5_files.join(' ')
    """
    mkdir -p "${outputdir}"
    python -m macproqc_helpers visualize -hdf5_files ${hdf5_files_list} -output "${outputdir}" -spikein_columns "${spikein_columns}" ${spike_ins_tab} $args

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
