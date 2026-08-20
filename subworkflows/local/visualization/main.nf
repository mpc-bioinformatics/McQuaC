include { QCVISUALIZATION } from '../../../modules/local/qcvisualization/main'

workflow VISUALIZATION {

    take:
    hdf5_files       // channel: [ val(meta), hdf5s]
    outputdir        // val: output directory
    spike_ins_table  // path: spike-ins table

    main:

    QCVISUALIZATION ( hdf5_files, outputdir, spike_ins_table )

    emit:
    jsons = QCVISUALIZATION.out.jsons           // channel: [ val(meta), [ json ] ]
    htmls = QCVISUALIZATION.out.htmls           // channel: [ val(meta), [ html ] ]
    output_tables = QCVISUALIZATION.out.output_tables // channel: [ val(meta), [ csv, tsv, xlsx ] ]
    ms1_maps = QCVISUALIZATION.out.ms1_maps     // channel: [ val(meta), ms1_maps ]
    additional_plots = QCVISUALIZATION.out.additional_plots // channel: [ val(meta), additional_plots ]
    bruker_calibrants = QCVISUALIZATION.out.bruker_calibrants // channel: [ val(meta), bruker_calibrants ]
}
