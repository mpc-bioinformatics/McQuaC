#!/usr/bin/env nextflow
nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'

workflow combine_metric_csvs {
	take:
		runbase_to_csvs
	
	main:
		headers_file = write_headers(runbase_to_csvs.first())
		datasection_files = write_data(runbase_to_csvs)
		
		final_table = combine_metric_csv_files(
			headers_file.concat(datasection_files).collect()
		)
	
	emit:
		final_table
}


process write_headers {
	container {python_image}

	input:
	tuple val(runBaseName), path(metrics)

	output:
	path "headers.csv"

	"""
	echo -n "file_and_analysis_timestamp" > headers.csv

	declare -a files=(${metrics})
	for file in \${files[@]}
	do
		# some files contain windows linebreaks, remove these
		line=\$(head -n -1 \${file} | sed -e "s;\\r\$;;g")
		echo -n ",\${line}" >> headers.csv
	done
	"""
}


process write_data {
	container {python_image}

	input:
	tuple val(runBaseName), path(metrics)

	output:
	path "${runBaseName}-data.csv"

	"""
	echo -n \$(date +"%Y-%m-%d_%H:%M:%S") > ${runBaseName}-data.csv

	declare -a files=(${metrics})
	for file in \${files[@]}
	do
		# some files contain windows linebreaks, remove these
		line=\$(tail -n 1 \${file} | sed -e "s;\\r\$;;g")
		echo -n ",\${line}" >> ${runBaseName}-data.csv
	done
	"""
}


process combine_metric_csv_files {
	container {python_image}

	input:
	path metrics_csvs

	output:
	path "quality_control.csv"

	"""
	declare -a files=(${metrics_csvs})
	for file in \${files[@]}
	do
		cat \$file >> quality_control.csv
		echo "" >> quality_control.csv
	done
	"""
}