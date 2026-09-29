/*
 * Create a ThermoRawFileParser (also alphatims-compatible) XIC extraction config
 * from a PSM mzTab file and a spike-ins CSV
 **/
process XICEXTRACTIONCONFIG {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/mpc-bioinformatics/macproqc-helpers:sha-7a956ec"

    input:
    tuple val(meta), path(psm_mztab_file)
    path spike_ins

    output:
    tuple val(meta), path("*.json"), path("*.ident.csv"), emit: config
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python -m macproqc_helpers create_thermorawfileparser_xic_config \\
        ${args} \\
        -icsv ${spike_ins} \\
        -iidents ${psm_mztab_file} \\
        -ojson ${psm_mztab_file}.json \\
        -oidentifications ${psm_mztab_file}.ident.csv

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

    echo '{}' > ${prefix}.json
    touch ${prefix}.ident.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
        macproqc-helpers: \$(python -m macproqc_helpers --version | sed 's;__main__.py ;;')
    END_VERSIONS
    """
}
