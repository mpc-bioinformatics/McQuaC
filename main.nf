#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* Mass centric Quality Control Workflow --> McQuaC

This is the QC-Workflow which generates various metrics from measured ThermoFischer-data ("RAW-files"). 
Here we extract XICs, Identifications and other various information and save them to a database.

The <TODO_VIS>.nf-Script then can be used to genereate various plots from the extracted data for inspection.

Example call:
nextflow run \
	-resume main.nf \
	--main_input_folder <Path_to_folder_of_raws>/raws \
	--main_fasta_file <Path_to_FASTA_file>.fasta 

*/

// Include all the needed workflows from the sub-workflows
// Extend this to also extend the QC-Workflow 
include {convert_raws_to_mzml} from './src/io/raw_file_conversion.nf'
include {identification_with_comet; identification_with_comet as identification_labelled_with_comet} from './src/identification/comet.nf'
include {pia_analysis_full; pia_analysis_psm_only; pia_extract_metrics} from './src/pia.nf'
include {retrieve_spike_ins_information} from './src/retrieve_spike_ins.nf'
include {get_feature_metrics} from './src/feature_detection.nf'
include {get_headers; get_mzml_infos} from './src/metrics/ms_run_metrics.nf'
include {combine_metric_hdf5} from './src/io/combine_metric_hdf5.nf'
include {output_processing_success} from './src/io/output_processing_success.nf'
include {visualization} from './src/visualization.nf'

// Parameters required for the standalone execution of this main-nextflow script
params.main_input_folder = "" // The folder containing the raw/.d files or hdf5 files (if visualize_only is true)
params.main_fasta_file = "" // A SINGLE-Fasta-file of the species to be searched (should also contain the SpikeIns if needed)
params.mcquac_params_file = "${baseDir}/example_configurations/mcquac_params.json" // the main parameters file for McQuaC

params.spike_ins_table = "${baseDir}/example_configurations/spike_ins.csv" // The information about spike-ins 
params.main_outdir = "./results"  // Output-Directory of the result files

// Parameter to only visualize the results (as hdf5 files) without running the whole workflow
params.visualize_only = false // If true, the workflow will only visualize the results (hdf5 files) in main_input_folder and not run the whole workflow.

// Parameters for visualization script
params.rt_unit = "sec" // Unit of the retention time, either sec for seconds or min for minutes.
params.output_column_order = "''" // Order of columns in the output table
params.spikein_columns = "Maximum_Intensity,RT_at_Maximum_Intensity,PSMs,Delta_to_expected_RT" // Columns of the spike-in dataframes that should end up in the result table
params.output_table_type = "csv" // Type of the output table, either csv or xlsx
params.height_barplots = 700 // Height of the barplots in pixels
params.width_barplots = 0 // Width of the barplots in pixels, 0 = flexible width
params.height_pca = 1000 // Height of the PCA plot in pixels
params.width_pca = 1000 // Width of the PCA plot in pixels
params.height_ionmaps = 10 // Height of the ion maps in inches
params.width_ionmaps = 10 // Width of the ion maps in inches

// Here are some optional Parameters which can be set if needed
params.search_spike_ins = true // Parameter to check if we execute a isa specific xic extraction (NOTE: FASTA has to contain the SpikeIns too!)
params.search_labelled_spikeins = true // Perform a special ID and look for labelled peptides

// Following are parameters to be passed to the sub-workflows

// Memory for the Thermo Raw File Parser, used 24 GB for a Raw file with 257409 MS scans 
/// and 4GB for a Raw file with 11352 MS scans (measured with `/usr/bin/time -v ...`). 10 GB seems legit for most cases.
params.ms_run_metrics__thermo_raw_mem = "10.GB"
// Memory for the tdf2mzml, used 0.39 GB for a Raw file with 298748 MS scans 
/// and 0.14GB for a Raw file with 35023 MS scans (measured with `/usr/bin/time -v ...`). 5 GB seems more then enough.
params.ms_run_metrics__bruker_raw_mem = "1.GB"
/// Tracing showed up to 4.7 GB virtual memory for 30000 MS scans
params.ms_run_metrics__mzml_mem = "10.GB"
// Set if you want to extract specific headers from Thermo raw measurements, otherwise the default is used.
// Have a look into the corresponding python script for the headers.
params.ms_run_metrics__thermo_headers = ""
 // Set if you want to extract specific headers from Bruker raw measurements, otherwise the default is used.
 // Have a look into the corresponding python script for the headers.
params.ms_run_metrics__bruker_headers = ""


// MAIN WORKFLOW
workflow {
	if (params.visualize_only) { 		 // only visualization
		list_hdf5_files = file(params.main_input_folder + '/*.hdf5')
		visualization(list_hdf5_files, params.main_outdir, params.rt_unit, params.output_column_order, params.spikein_columns, params.output_table_type, params.search_spike_ins)
  }
	else { // whole workflow execution

    // Retrieve input files
    thermo_raw_files = Channel.fromPath(params.main_raw_spectra_folder + "/*.raw")
    bruker_raw_folders = Channel.fromPath(params.main_raw_spectra_folder + "/*.d", type: 'dir')

    main_outdir = Channel.fromPath(params.main_outdir).first()

    // .first() convert the queue channel with only one file to a value channel, making it possible to use multiple time
    // e.g. to automatically start multiple concurrent identifications (no need for map each raw file with the fasta and config file)
    fasta_file = Channel.fromPath(params.main_fasta_file).first()
    mcquac_params_file = Channel.fromPath(params.mcquac_params_file).first()

    raw_files = thermo_raw_files.concat(bruker_raw_folders)
    // File conversion into open formats
    mzmls = convert_raws_to_mzml(thermo_raw_files, bruker_raw_folders)

    // Retreive MZML Metrics
    mzml_metrics = get_mzml_infos(mzmls, mcquac_params_file, params.ms_run_metrics__mzml_mem)

    // Identify spectra using Comet
    comet_ids = identification_with_comet(mzmls, fasta_file, mcquac_params_file, false, main_outdir)

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
      comet_labelled_ids = identification_labelled_with_comet(mzmls, fasta_file, mcquac_params_file, true, main_outdir)

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
		feature_metrics = get_feature_metrics(mzmls, pia_report_psm_mztabs, mcquac_params_file)

	  // Get Thermo/Bruker specific information from raw_spectra
	  custom_header_infos = get_headers(thermo_raw_files, params.ms_run_metrics__thermo_raw_mem, params.ms_run_metrics__thermo_headers,
		  	bruker_raw_folders, params.ms_run_metrics__bruker_raw_mem, params.ms_run_metrics__bruker_headers)

		// Concatenate to one merged metric CSV
		hdf5s_per_run = mzml_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-mzml_info.hdf5')), file)}
		if (params.search_spike_ins) {
			hdf5s_per_run = hdf5s_per_run.concat(
				spike_in_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-spikeins.hdf5')), file)}
			)
		}
		hdf5s_per_run = hdf5s_per_run
			.concat(feature_metrics.map{file -> tuple(file.name.take(file.name.lastIndexOf('-features.hdf5')), file)})
			.concat(pia_extract_csv.map{file -> tuple(file.name.take(file.name.lastIndexOf('-pia_extraction.hdf5')), file)})
			.concat(custom_header_infos.map{file -> tuple(file.name.take(file.name.lastIndexOf('-custom_headers.hdf5')), file)})
			.groupTuple()

    combined_metrics = combine_metric_hdf5(hdf5s_per_run, main_outdir)

		// Visualize the results (and move them to the results folder)
		visualization(combined_metrics, params.main_outdir, params.rt_unit, params.output_column_order, params.spikein_columns, params.output_table_type, params.search_spike_ins, 
			params.height_barplots, params.width_barplots, params.height_pca, params.width_pca, params.height_ionmaps, params.width_ionmaps)

		output_processing_success(raw_files, hdf5s_per_run.toList().transpose().first().flatten())

	}
}

