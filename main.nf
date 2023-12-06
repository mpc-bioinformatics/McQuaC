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
PROJECT_DIR = workflow.projectDir
include {convert_to_mzml; convert_to_mgf; convert_to_idxml} from PROJECT_DIR + '/file_conversion.nf'
include {get_various_mzml_infos} from PROJECT_DIR + '/get_mzml_chromatogram_and_more.nf'
include {ident_via_comet} from PROJECT_DIR + '/identification_via_comet.nf'
include {execute_pia} from PROJECT_DIR + '/pia.nf'
include {retrieve_spikeins} from PROJECT_DIR + '/retrieve_spike_ins.nf' // We could also consider to expose params.spk_spike_ins, however it is always fixed for our ISA-stadard!
include {get_features} from PROJECT_DIR + '/get_features_in_raws.nf'
include {get_custom_headers} from PROJECT_DIR + '/get_custom_columns_from_file_directly.nf'
// Each script has its own UNIQUE-param-attribute and can be fine-tuned from this main.nf-script.
// The requiered params are also exposed in this script and are listed below:

// Parameters required for the standalone execution of this main-nextflow script
params.main_raw_spectra_folder = "" // The folder containing the raw spectra
params.main_fasta_file = "" // A SINGLE-Fasta-file of the species to be searched (should also contain the SpikeIns if needed)
params.main_comet_params = "${baseDir}/example_configurations/comet_config.txt" // Main-Search-Parameters for the comet search engine
params.main_outdir = "$PWD/results"  // Output-Directory of the Identification Results. Here it is <Input_File>.mzid


// Here are some optional Parameters which can be set if needed
params.main_is_isa = true // Parameter to check if we execute a isa specific xic extraction (NOTE: FASTA has to contain the SpikeIns too!)

// MAIN WORKFLOW
workflow {
	// Retrieve inpit files
	thermo_raw_files = Channel.fromPath(params.main_raw_spectra_folder + "/*.raw")
	bruker_raw_files = Channel.fromPath(params.main_raw_spectra_folder + "/*.d", type: 'dir')
	// .first() convert the queue channel with only one file to a value channel, making it possible to use multiple time
	// e.g. to automatically start multiple concurrent identifications (no need for map each raw file with the fasta and config file)
	fasta_file = Channel.fromPath(params.main_fasta_file).first()
	comet_params = Channel.fromPath(params.main_comet_params).first()


	// File conversion into open formats
	mzmls = convert_to_mzml(thermo_raw_files, bruker_raw_files)
	// mgfs = convert_to_mgf(mzmls)
	
	// // Retreive MZML Statistics
	get_various_mzml_infos(mzmls)

	// /* Identify with multiple search engines */
	// // Comet
	pepxmls = ident_via_comet(mzmls, fasta_file, comet_params)
	idxmls = convert_to_idxml(pepxmls)
	
	// // MS-GF+
	// // convert_to_mgf(rawspectra) //MS-GF+ Workflow requires MGFs
	// // TODO ||| HERE goes the code!
	// /* END Identify with multiple search engines */

	// // Execute PIA and filter by FDR ||| TODO do we combine the search engine results?
	execute_pia(idxmls)

	// Specific to ISA: Do XIC-Extraction if specified
	if (params.main_is_isa) {
		raw_files = thermo_raw_files.concat(bruker_raw_files)
		retrieve_spikeins(raw_files, execute_pia.out[0].map { it[0] })
	}

	// Run Feature Finding and Statistics
	get_features(mzmls, execute_pia.out[0].map { it[0] })

	// Get Thermospecific information from raw
	get_custom_headers(thermo_raw_files, bruker_raw_files)


	// // Concatenate to large csv
	combined_csvs = get_various_mzml_infos.out.collect().concat(
		retrieve_spikeins.out.collect(),
		get_features.out.map { it[1] }.collect(),
		execute_pia.out[1].collect(),
		get_custom_headers.out.collect()
	).collect()
	combine_output_to_table(combined_csvs)
	

	// // Visualize the results
	visualize_results(combine_output_to_table.out)

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
	QC_visualization.py -csv_file $complete_csv -output "." 
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

