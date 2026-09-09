include { QCVISUALIZATION } from '../../../modules/local/qcvisualization/main'

workflow VISUALIZATION {

    take:
    hdf5_files       // channel: [ val(meta), hdf5s]
    spike_ins_table  // path: spike-ins table

    main:

    def spikein_columns = params.spikein_columns ?: 'MS1 feature maximum intensity,retention time,count of identified spectra,Delta_to_expected_RT'

    QCVISUALIZATION (
        hdf5_files,
        spike_ins_table,
        params.figure_format,
        params.output_table_type,
        params.spikeins,
        params.rt_unit,
        params.output_column_order ?: '',
        spikein_columns,
        params.height_barplots,
        params.width_barplots,
        params.height_pca,
        params.width_pca,
        params.height_ionmaps,
        params.width_ionmaps,
    )

    emit:
    jsons = QCVISUALIZATION.out.jsons           // channel: [ val(meta), [ json ] ]
    htmls = QCVISUALIZATION.out.htmls           // channel: [ val(meta), [ html ] ]
    output_tables = QCVISUALIZATION.out.output_tables // channel: [ val(meta), [ csv, tsv, xlsx ] ]
    ms1_maps = QCVISUALIZATION.out.ms1_maps     // channel: [ val(meta), ms1_maps ]
    additional_plots = QCVISUALIZATION.out.additional_plots // channel: [ val(meta), additional_plots ]
    bruker_calibrants = QCVISUALIZATION.out.bruker_calibrants // channel: [ val(meta), bruker_calibrants ]
}
