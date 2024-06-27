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
include {identification_with_comet} from workflow.projectDir + '/src/identification/comet.nf'
include {identification_with_comet as identification_labelled_with_comet} from workflow.projectDir + '/src/identification/comet.nf'
include {convert_to_idxml} from workflow.projectDir + '/file_conversion.nf'
include {convert_to_idxml as convert_to_labelled_idxml} from workflow.projectDir + '/file_conversion.nf'
include {get_various_mzml_infos} from workflow.projectDir + '/get_mzml_chromatogram_and_more.nf'
include {pia_analysis_full; pia_analysis_psm_only; pia_extract_metrics} from workflow.projectDir + '/src/pia.nf'
include {retrieve_spike_ins_information} from workflow.projectDir + '/src/retrieve_spike_ins.nf'
include {get_features} from workflow.projectDir + '/get_features_in_raws.nf'
include {get_custom_headers} from workflow.projectDir + '/get_custom_columns_from_file_directly.nf'

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

	feature_finder_params = get_feature_finder_params_from_comet_params(comet_params)

	raw_files = thermo_raw_files.concat(bruker_raw_folders)
	// File conversion into open formats
	mzmls = convert_raws_to_mzml(thermo_raw_files, bruker_raw_folders)
	
	// Retreive MZML Statistics
	get_various_mzml_infos(mzmls)

	// Identify spectra using Comet
	comet_pepxmls = identification_with_comet(mzmls, fasta_file, comet_params, false)
	comet_idxmls = convert_to_idxml(comet_pepxmls)

	// Execute protein inference and filter by FDR
	pia_report_files = pia_analysis_full(comet_idxmls)
	pia_extract_csv = pia_extract_metrics(pia_report_files)

	// search additionally for labelled PSMs
	if (params.search_labelled_spikeins) {
		comet_labelled_pepxmls = identification_labelled_with_comet(mzmls, fasta_file, comet_params, true)
		comet_labelled_idxmls = convert_to_labelled_idxml(comet_labelled_pepxmls)

		// set the filter to true to count only FDR filtered labelled PSMs - but the FDR is skewed anyways, as the label is set to static!
		labelled_pia_report_files = pia_analysis_psm_only(comet_labelled_idxmls, false)
	}

	// extract spike-ins information
	if (params.search_spike_ins) {
		spike_ins_table = Channel.fromPath(params.spike_ins_table).first()

		if (params.search_labelled_spikeins) {
			psm_results = labelled_pia_report_files
		} else {
			psm_results = pia_report_files
				.toList()
            	.transpose()
            	.first()
            	.flatten()
		}

		retrieve_spike_ins_information(raw_files, psm_results, spike_ins_table)
	}
	 
	// // Run Feature Finding and Statistics
	// get_features(mzmls, pia_analysis.out[0].map { it[0] }, feature_finder_params)

	// // Get Thermospecific information from raw
	// get_custom_headers(thermo_raw_files, bruker_raw_folders)


	// // // Concatenate to large csv
	// combined_csvs = get_various_mzml_infos.out.collect()
	// if (params.main_is_isa) {
	// 	combined_csvs = combined_csvs.concat(
	// 		retrieve_spikeins.out.collect(),
	// 		get_features.out.map { it[1] }.collect(),
	// 		execute_pia.out[1].collect(),
	// 		get_custom_headers.out.collect()
	// 	).collect()
	// } else {
	// 	combined_csvs = combined_csvs.concat(
	// 		get_features.out.map { it[1] }.collect(),
	// 		execute_pia.out[1].collect(),
	// 		get_custom_headers.out.collect()
	// 	).collect()
	// }
	// combine_output_to_table(combined_csvs)
	

	// // // Visualize the results
	// visualize_results(combine_output_to_table.out)

}

/**
 * This process is used to convert the comet parameters into the parameters needed for the feature finder.
 * 
 * @param comet_params The comet parameters file
 * @return The feature finder parameters (value channel)
 */
process get_feature_finder_params_from_comet_params {
	container 'mpc/nextqcflow-python:latest'

	input:
	path comet_params

	output:
	stdout

	"""
	comet_params_to_feature_finder_params.py -c ${comet_params}
	"""
}

process combine_output_to_table {
	container 'mpc/nextqcflow-python:latest'

	publishDir "${params.main_outdir}/qc_results", mode:'copy'

	input:
	path(input_csv_files)

    output:
    path("quality_control.csv")

    """
	CONCAT_CSVS=""
    for file in $input_csv_files
    do
        CONCAT_CSVS+="\$file,"
    done
    CONCAT_CSVS=\$(echo \$CONCAT_CSVS | rev | cut -c2- | rev)

	unify_csv_tables.py -out_csv quality_control.csv -input_csvs \$CONCAT_CSVS
    """
}

process visualize_results {
	container 'mpc/nextqcflow-python:latest'

	publishDir "${params.main_outdir}/qc_results", mode:'copy'

	input:
	path(complete_csv)

    output:
    path("*.json")
	path("*.html")
	path("*.csv")
	path("fig13_ionmaps")
	path("THERMO_PLOTS_FIG15")

	"""
	if ${params.main_is_isa}
	then 
		QC_visualization.py -csv_file $complete_csv -output "." -isa
	else
		QC_visualization.py -csv_file $complete_csv -output "."
	fi
    """
}

// // TODO Usefull bits: 
// params.help = false

// if(params.help) {
// 	println(""" \nusage : ~/nextflow main.nf [operators] --input inputFiles
// 		example : ~/nextflow main.nf --input \"path/fastas/\"
// 		operators: --help calls this help option\n """)
// 	System.exit(0)
// 	} else {
// 		dir= workflow.projectDir.getParent() + "/results/"

// include {msgfplus_buildsa; msgfplus_search_mgf} from PROJECT_DIR + "/identification_via_msgfplus.nf"
// workflow msgfplus{
//     // Get all MGF files which should be identified
//  //   mgfs = Channel.fromPath(converter.out)

//     // Get FASTA-file
//     fasta_file = Channel.fromPath(params.fasta_file)

//     // Get Modification Parameters file
//     modifications_file = Channel.fromPath(params.search_parameter_file)

// 	take: data

// 	main: 
//     // Build indexed fasta for MSGFPLUS
//     	msgfplus_buildsa(fasta_file)
//     // Combined channel replicated the indexed fasta for each MGF to be reused
//     	combined_channel = fasta_file
//         	.combine(modifications_file)
//         	.combine(data)
//         	.combine(msgfplus_buildsa.out.toList())
// 	// Start search
//     	msgfplus_search_mgf(combined_channel)
// 		msgfplus_search_mgf.out.view()

// 	emit:
// 		msgfplus_search_mgf.out
// }

