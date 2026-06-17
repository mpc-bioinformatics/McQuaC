/*
 * Convert Bruker .d-Folder to mzML using tdf2mzml
 **/
process TDF2MZML {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE'
        : 'quay.io/medbioinf/tdf2mzml:0.4'}"

    input:
    tuple val(meta), path(d_folder)

    output:
    tuple val(meta), path("*.mzML"), emit: spectra
    tuple val("${task.process}"), val('tdf2mzml'), val("0.4"), topic: versions, emit: versions_tdf2mzml

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    export MKL_NUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    tdf2mzml ${args} -i ${d_folder} --compression "zlib" -o ${prefix}.mzML
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.mzML
    """
}
