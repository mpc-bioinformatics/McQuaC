#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* MPC-Nextflow-Quality-Control-Workflow --> MPC-QuaC-Workflow

This is the QC-Workflow which generates various statistics from measured ThermoFischer-data ("RAW-files"). 
Here we extract XICs, Identifications and other various information and save them to a database.

The <TODO_VIS>.nf-Script then can be used to genereate various plots from the extracted data for inspection.

Example call:
nextflow run \
	-resume main.nf \
	--main_raw_spectra_folder <Path_to_folder_of_raws>/raws \
	--main_fasta_file <Path_to_FASTA_file>.fasta 

*/

// Include all the needed workflows from the sub-workflows
// Extend this to also extend the QC-Workflow 
include {convert_raws_to_mzml} from workflow.projectDir + '/src/io/raw_file_conversion.nf'
include {identification_with_comet; identification_with_comet as identification_labelled_with_comet} from workflow.projectDir + '/src/identification/comet.nf'
include {pia_analysis_full; pia_analysis_psm_only; pia_extract_metrics} from workflow.projectDir + '/src/pia.nf'
include {retrieve_spike_ins_information} from workflow.projectDir + '/src/retrieve_spike_ins.nf'
include {get_feature_metrics} from workflow.projectDir + '/src/feature_detection.nf'
include {get_headers; get_mzml_infos} from workflow.projectDir + '/src/metrics/ms_run_metrics.nf'
include {combine_metric_csvs} from workflow.projectDir + '/src/io/combine_metric_csvs.nf'
include {output_processing_success} from workflow.projectDir + '/src/io/output_processing_success.nf'

// Parameters required for the standalone execution of this main-nextflow script
params.main_raw_spectra_folder = "" // The folder containing the raw spectra
params.main_fasta_file = "" // A SINGLE-Fasta-file of the species to be searched (should also contain the SpikeIns if needed)
params.main_comet_params = "${baseDir}/example_configurations/high-high.comet.params" // Main-Search-Parameters for the comet search engine
params.spike_ins_table = "${baseDir}/example_configurations/spike_ins.csv" // The information about spike-ins 
params.main_outdir = "$PWD/results"  // Output-Directory of the Identification Results. Here it is <Input_File>.mzid

// Here are some optional Parameters which can be set if needed
params.search_spike_ins = true // Parameter to check if we execute a isa specific xic extraction (NOTE: FASTA has to contain the SpikeIns too!)
params.search_labelled_spikeins = true // Perform a special ID and look for labelled peptides


// MAIN WORKFLOW
workflow {
	// Retrieve input files
	thermo_raw_files = Channel.fromPath(params.main_raw_spectra_folder + "/*.raw")
	bruker_raw_folders = Channel.fromPath(params.main_raw_spectra_folder + "/*.d", type: 'dir')

	// .first() convert the queue channel with only one file to a value channel, making it possible to use multiple time
	// e.g. to automatically start multiple concurrent identifications (no need for map each raw file with the fasta and config file)
	fasta_file = Channel.fromPath(params.main_fasta_file).first()
	comet_params = Channel.fromPath(params.main_comet_params).first()

	raw_files = thermo_raw_files.concat(bruker_raw_folders)
	// File conversion into open formats
	mzmls = convert_raws_to_mzml(thermo_raw_files, bruker_raw_folders)
	
	// Retreive MZML Metrics
	mzml_metrics = get_mzml_infos(mzmls)

	// Identify spectra using Comet
	comet_ids = identification_with_comet(mzmls, fasta_file, comet_params, false)

	// Execute protein inference and filter by FDR
	pia_report_files = pia_analysis_full(comet_ids.mzids)
	pia_report_psm_mztabs = pia_report_files
				.toList()
            	.transpose()
            	.first()
            	.flatten()
	pia_extract_csv = pia_extract_metrics(pia_report_files)

	// search additionally for labelled PSMs
	if (params.search_labelled_spikeins) {
		comet_labelled_ids = identification_labelled_with_comet(mzmls, fasta_file, comet_params, true)

		// set the filter to true to count only FDR filtered labelled PSMs - but the FDR is skewed anyways, as the labelling is set to "static"!
		labelled_pia_report_files = pia_analysis_psm_only(comet_labelled_ids.mzids, false)
	}

	// extract spike-ins information
	if (params.search_spike_ins) {
		spike_ins_table = Channel.fromPath(params.spike_ins_table).first()

		if (params.search_labelled_spikeins) {
			psm_results = labelled_pia_report_files
		} else {
			psm_results = pia_report_psm_mztabs
		}

		spike_in_metrics = retrieve_spike_ins_information(raw_files, psm_results, spike_ins_table)
	}
	 
	// Run Feature Finding
	feature_metrics = get_feature_metrics(mzmls, pia_report_psm_mztabs, comet_params)

	// Get Thermospecific information from raw
	custom_header_infos = get_headers(thermo_raw_files, bruker_raw_folders)

	// Concatenate to one merged metric CSV
	csvs_per_run = mzml_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-mzml_info.csv')), file)}
	if (params.search_spike_ins) {
		csvs_per_run = csvs_per_run.concat(
			spike_in_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-spikeins.csv')), file)}
		)
	}
	csvs_per_run = csvs_per_run
		.concat(feature_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-features.csv')), file)})
		.concat(pia_extract_csv.map{file -> tuple(file.name.take(file.name.lastIndexOf('-pia_extraction.csv')), file)})
		.concat(custom_header_infos.map{file -> tuple(file.name.take(file.name.lastIndexOf('-custom_headers.csv')), file)})
		.groupTuple()
	
	combined_metrics = combine_metric_csvs(csvs_per_run)

	// Visualize the results (and move them to the results folder)
	visualize_results(combined_metrics)

	output_processing_success(raw_files, csvs_per_run.toList().transpose().first().flatten())
}

/**
 * This process calls the script to finally draw the visualizations
 * 
 * @param comet_params The comet parameters file
 * @return All the visualizations (HTMLs, JSONs, CSVs)
 */
process visualize_results {
	container 'mpc/nextqcflow-python:latest'

	publishDir "${params.main_outdir}/qc_results", mode:'copy'		// TODO: this should probably rather use the new reporting facilities

	input:
	path(combined_metrics_csv)

    output:
    path("*.json")
	path("*.html")
	path("*.csv")
	path("fig13_ionmaps")
	path("THERMO_PLOTS_FIG15")

	"""
	if ${params.search_spike_ins}
	then 
		QC_visualization.py -csv_file ${combined_metrics_csv} -output "." -spikeins
	else
		QC_visualization.py -csv_file ${combined_metrics_csv} -output "."
	fi
    """
}
