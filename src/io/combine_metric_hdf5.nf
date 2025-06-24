#!/usr/bin/env nextflow

/**
 * Workflows and processes for merging the generated quality control metrics
 */


nextflow.enable.dsl=2

/*
 * Combines the metrics into one big HDF5 file
 * 
 * @param runbase_to_hdf5 tuple of val(fileBaseName) and list(path(metric_HDF5s))
 *
 * @return final_table a merged, big HDF5 file containing data of all files
 */
workflow combine_metric_hdf5 {
	take:
		runbase_to_hdf5
	
	main:
	 	merged_files = merge_metrics(runbase_to_hdf5)
		complete_hdf_files = merged_files.collect()
	
	emit:
		complete_hdf_files
}

/*
 * Writes the metrics for one raw file into one hdf5 file
 * 
 * @param runbase_to_hdf5 tuple of val(fileBaseName) and list(path(metric_hdf5s))
 *
 * @return final_table a merged, big hdf5 file of the single metrics
 */

process merge_metrics {
	container {python_image}

	publishDir "${params.main_outdir}/qc_hdf5_data", mode:'copy'		// TODO: this should probably rather use the new reporting facilities

	input:
	tuple val(runBaseName), path(metrics)

	output:
	path "${runBaseName}.hdf5"

	"""
	combine_hdf5_files.py -hdf_out_name ${runBaseName}.hdf5 $metrics
	"""
}
