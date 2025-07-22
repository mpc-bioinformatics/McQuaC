#!/usr/bin/env nextflow

/**
 * Workflows for vizualization of the reults
 */

nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'

workflow visualization {
    take: 
        combined_metrics
        main_outdir
        rt_unit
        output_column_order
        spikein_columns
        output_table_type
        search_spike_ins

    main:
        visualize_results(combined_metrics, main_outdir, rt_unit, output_column_order, spikein_columns, output_table_type, search_spike_ins)

    emit:
        visualize_results.out[0]
        visualize_results.out[1]
        visualize_results.out[2]
        visualize_results.out[3]
        visualize_results.out[5]
        visualize_results.out[6]
}

process visualize_results {
    container { python_image }

    publishDir "${main_outdir}/qc_results", mode:'copy'		// TODO: this should probably rather use the new reporting facilities

    input:
    path combined_metrics
    val main_outdir
    val rt_unit
    val output_column_order
    val spikein_columns
    val output_table_type
    val search_spike_ins
    


    output:
    path("*.json")
    path("*.html")
    path("*.${output_table_type}")
    path(combined_metrics)
    path("fig13_MS1_map")
    path("fig15_additional_headers")
    path("fig16_BRUKER_calibrants")

    """
    if ${search_spike_ins}
    then 
      QC_visualization.py -hdf5_files ${combined_metrics} -output "." -spikeins -RT_unit ${rt_unit} -output_column_order ${output_column_order} -spikein_columns ${spikein_columns} -output_table_type ${output_table_type}
    else
      QC_visualization.py -hdf5_files ${combined_metrics} -output "." -RT_unit ${rt_unit} -output_column_order ${output_column_order} -spikein_columns ${spikein_columns} -output_table_type ${output_table_type}
    fi
      """
}
