include { OPENMS_DECOYDATABASE } from '../../../modules/nf-core/openms/decoydatabase/main'
include { COMETCONFIG; COMETCONFIG as COMETCONFIG_LABELLED } from '../../../modules/local/cometconfig/main'

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

    // create fiel channel for comet config template, either user-specified or default
    ch_comet_config_template = comet_config_template
        ? channel.fromPath(comet_config_template, checkIfExists: true)
        : channel.fromPath("${projectDir}/assets/default_configs/comet.params", checkIfExists: true)

    // add an id to the comet config template channel
    ch_comet_config_template = ch_comet_config_template.map { params ->
        def params_id = comet_config_template
            ? 'user-template'
            : 'default'
        [ [id: params_id], params]
    }

    // create the config file (for unlabelled search)
    COMETCONFIG(
        ch_comet_config_template
    )

    // combine the channels for unlabelled search
    ch_comet_input = ch_mzml.combine(ch_fasta).combine(COMETCONFIG.out.params)
        .map { mzml_meta, mzml, fasta_meta, fastafile, params_meta, comet_params ->
            def meta = mzml_meta + [mzml_id: mzml_meta.id, label_search: false, fasta_id: fasta_meta.id, params_id: params_meta.id]
            meta.id = "${meta.mzml_id}-unlabelled"
            [meta, mzml, fastafile, comet_params]
        }

    // TODO: change to if (search_label == true)....
    if (true) {
        // create the config file (for labelled search)
        COMETCONFIG_LABELLED(
            ch_comet_config_template
        )

        // create channels for labelled search
        ch_comet_input_labelled = ch_mzml.combine(ch_fasta).combine(COMETCONFIG_LABELLED.out.params)
            .map { mzml_meta, mzml, fasta_meta, fastafile, params_meta, comet_params ->
                def meta = mzml_meta + [mzml_id: mzml_meta.id, label_search: true, fasta_id: fasta_meta.id, params_id: params_meta.id]
                meta.id = "${meta.mzml_id}-labelled"
                [meta, mzml, fastafile, comet_params]
            }

        // mix all params for comet call
        ch_comet_input = ch_comet_input.mix(ch_comet_input_labelled)
    }

    // emit:
    // emit outputs, as soon as created
}
