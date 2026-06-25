include { OPENMS_DECOYDATABASE } from '../../../modules/nf-core/openms/decoydatabase/main'

workflow IDENT_DDA {

    take:
    fasta                   // val: fasta path
    skip_decoy_generation   // val: true to skip the decoy generation
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

    // comet id
    // needs: params-file (global for now)
    // get inputs like before...

    //emit:

}
