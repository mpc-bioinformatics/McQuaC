
process QCVISUALIZATION {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'ghcr.io/mpc-bioinformatics/macproqc-helpers:latest'}"

    input:
    tuple val(meta), path(hdf5_files)
    val outputdir
    path spike_ins_table

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("${outputdir}/*.json"), emit: jsons, optional: true
    tuple val(meta), path("${outputdir}/*.html"), emit: htmls, optional: true
    tuple val(meta), path("${outputdir}/*.{csv,tsv,xlsx}"), emit: output_tables, optional: true
    tuple val(meta), path("${outputdir}/fig13_MS1_map"), emit: ms1_maps, optional: true
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
    // --meta "${meta}"

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/624977dfaf562211e68a8a868ca80acc8461f1ac/modules/nf-core/cutadapt/main.nf#L34-L46
    //               Complex example: https://github.com/nf-core/modules/blob/88d43dad73a675e66bff49ebb57fe657a5909018/modules/nf-core/bedtools/split/main.nf#L32-L43
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args
    
    touch ${prefix}.jsons

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
