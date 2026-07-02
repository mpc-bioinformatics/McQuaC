include { OPENMS_DECOYDATABASE } from '../../../modules/nf-core/openms/decoydatabase/main'
include { COMETCONFIG } from '../../../modules/local/cometconfig/main'

workflow IDENT_DDA {
    take:
    fasta                   // val: fasta path
    skip_decoy_generation   // val: true to skip the decoy generation
    comet_config_template   // val: path to comet config template (or null to use default)
    ch_mzml                 // channel: [ val(meta), mzmls ]

    main:

    // for now, only one FASTA is used for all runs
    ch_fasta = channel.fromPath(fasta, checkIfExists: true)
        .map { f -> [[id: "global_fasta"], f] }

    if (!skip_decoy_generation) {
        OPENMS_DECOYDATABASE(
            ch_fasta
        )
        ch_fasta = OPENMS_DECOYDATABASE.out.decoy_fasta
    }

    ch_comet_config_template = comet_config_template
        ? channel.fromPath(comet_config_template, checkIfExists: true)
        : channel.fromPath("${projectDir}/assets/default_configs/comet.params", checkIfExists: true)

    ch_comet_config_template = ch_comet_config_template.map { params ->
        def params_id = comet_config_template
            ? 'user-template'
            : 'default'
        [ [id: params_id], params]
    }

    COMETCONFIG(
        ch_comet_config_template
    )
}
