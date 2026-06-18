// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { QCVISUALIZATION } from '../../../modules/local/qcvisualization/main'

workflow VISUALIZATION {

    take:
    // TODO nf-core: edit input (take) channels
    hdf5_files // channel: [ val(meta), [ bam ] ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    QCVISUALIZATION ( hdf5_files )

    emit:
    // TODO nf-core: edit emitted channels
    jsons = QCVISUALIZATION.out.jsons           // channel: [ val(meta), [ bam ] ]
}
