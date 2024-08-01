#!/usr/bin/env nextflow

/**
 * Workflows and processes for merging the generated quality control metrics
 */


nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'

/*
 * Combines the metrics into one big CSV file
 * 
 * @param runbase_to_csvs tuple of val(fileBaseName) and list(path(metric_CSVs))
 *
 * @return final_table a merged, big CSV file containing data of all files
 */
workflow combine_metric_csvs {
	take:
		runbase_to_csvs
	
	main:
	 	merged_files = merge_metrics(runbase_to_csvs)
		final_table = merged_files.collectFile(name: "quality_control.csv", keepHeader: true)
	
	emit:
		final_table
}

/*
 * Writes the metrics for one raw file into one CSV file (with one line)
 * 
 * @param runbase_to_csvs tuple of val(fileBaseName) and list(path(metric_CSVs))
 *
 * @return final_table a merged, big CSV file of teh single metrics
 */

process merge_metrics {
	container {python_image}

	input:
	tuple val(runBaseName), path(metrics)

	output:
	path "${runBaseName}-data.csv"

	"""
	csvfile="${runBaseName}-data.csv"

	# first write the headers
	echo -n "file_and_analysis_timestamp" > \${csvfile}

	declare -a files=(${metrics})
	for file in \${files[@]}
	do
		# some files contain windows linebreaks, remove these
		line=\$(head -n -1 \${file} | sed -e "s;\\r\$;;g")
		echo -n ",\${line}" >> \${csvfile}
	done

	# next line, write out the data
	echo >> \${csvfile}

	# first write out the current file_and_analysis_timestamp
	echo -n "${runBaseName}_" >> \${csvfile}
	echo -n \$(date +"%Y_%m_%d_%H_%M_%S") >> \${csvfile}

	declare -a files=(${metrics})
	for file in \${files[@]}
	do
		# some files contain windows linebreaks, remove these
		line=\$(tail -n 1 \${file} | sed -e "s;\\r\$;;g")
		echo -n ",\${line}" >>  \${csvfile}
	done

	# final line break (makes things easier)
	echo >> \${csvfile}
	"""
}