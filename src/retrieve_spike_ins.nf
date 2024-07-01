#!/usr/bin/env nextflow
nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'
thermorawfileparser_image = 'quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0'

params.max_parallel_xic_extractors = Runtime.runtime.availableProcessors() / 2

/*
 * Extracts for each defined spike-in the XIC in the respective m/z and RT region, looks
 * for identifications of the respective sequence and returns the highest peak of teh XIC,
 * its RT and the number of identifications.
 *
 * @param raw_files the raw files (.raw or .d)
 * @param psm_mztab_files the PSM results in mzTab
 * @param spike_ins_table the spike-ins defined in CSV format
 *
 * @return spike_in_metrics the metrics of teh spike-ins
 */
workflow retrieve_spike_ins_information {
    take:
        raw_files
        psm_mztab_files
        spike_ins_table
    
    main:
        // Finally, we generate the input json, retrieve it via trfp and parse back this results into a csv-format
        json_and_identifications = generate_json_and_identifications(psm_mztab_files, spike_ins_table)

        // Map the run basename to the raw files and to their corresponding mzTab files
        runbase_to_raw_and_json = 
            raw_files.map {
                file -> tuple(file.baseName, file)
            }.join(
                json_and_identifications.json.map{
                    file -> tuple(file.name.take(file.name.lastIndexOf('-trfp_input.json')), file)
                },
                by: 0
            )
        
        // branch to Thermo's .raw and Bruker's .d
        runbase_to_raw_and_json.branch {
            thermo: it[1].getExtension() == 'raw'
            bruker: it[1].getExtension() == 'd'
        }.set{ branched_runbase_to_raw_and_json }

        thermo_xics = retrieve_xics_from_thermo_raw_spectra(branched_runbase_to_raw_and_json.thermo)
        bruker_xics = retrieve_xics_from_bruker_raw_spectra(branched_runbase_to_raw_and_json.bruker)
        
        xics = thermo_xics.concat(bruker_xics)

        runbase_to_xics_and_identifications = 
            xics.map {
                file -> tuple(file.baseName, file)
            }.join(
                json_and_identifications.identifications.map{
                    file -> tuple(file.name.take(file.name.lastIndexOf('-identifications.csv')), file)
                },
                by: 0
            )

        spike_in_metrics = get_spike_in_metrics(runbase_to_xics_and_identifications, spike_ins_table)

    emit:
        spike_in_metrics
}

/*
 * Creates a json file for the XIC extraction of the spike-ins. For this, the
 * spike-ins file and the identification results are used to refine the retention
 * times for the extraction.
 *
 * @param psm_mztab_files the identifications in the PIA generated mzTab
 * @param spike_ins the spike-ins file
 *
 * @return *-trfp_input.json the JSON for the XIC extraction
 * @return *-identifications.csv a file mapping from the sequences to found identifications
 */
process generate_json_and_identifications {
    container { python_image }

    input:
    path psm_mztab_files
    path spike_ins

    output:
    path "${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-trfp_input.json", emit: json
    path "${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-identifications.csv", emit: identifications

    """
    create_spike_in_xic_json.py -icsv ${spike_ins} -iidents ${psm_mztab_files} -ojson ${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-trfp_input.json -oidentifications ${psm_mztab_files.name.take(psm_mztab_files.name.lastIndexOf('-piaExport-PSM.mzTab'))}-identifications.csv
    """
}

/**
 * Performs the XIC extraction on Thermo raw files
 *
 * @param run_basename (not used)
 * @param raw Thermo raw file
 * @param xic_json the JSON for the XIC extraction
 *
 * @return baseName.json the extracted XICs in JSON format
 */ 
process retrieve_xics_from_thermo_raw_spectra {
    container { thermorawfileparser_image }

    maxForks params.max_parallel_xic_extractors

    input:
    tuple val(run_basename), path(raw), path(xic_json)

    output:
    path "${raw.baseName}.json"

    """
    thermorawfileparser xic -i ${raw} -j ${xic_json}
    """
}

/**
 * Performs the XIC extraction on Bruker .d files
 *
 * @param run_basename (not used)
 * @param raw Bruker .d file
 * @param xic_json the JSON for the XIC extraction
 *
 * @return baseName.json the extracted XICs in JSON format
 */ 
process retrieve_xics_from_bruker_raw_spectra {
    container { python_image }

    maxForks params.max_parallel_xic_extractors

    input:
    tuple val(run_basename), path(raw), path(xic_json)

    output:
    path "${raw.baseName}.json"

    """
    extract_xic_bruker.py -d_folder ${raw} -in_json ${xic_json} -out_json ${raw.baseName}.json
    """
}

/**
 * Creates the metrics of the spike in extraction, regarding XICs, identification and RT deltas
 *
 * @param run_basename (not used)
 * @param xic_json the extracted XICs in JSON
 * @param identifications a CSV mapping the sequences to number of identifications
 *
 * @return baseName-spikeins.csv the metrics of the spike-ins extraction, in CSV
 */ 
process get_spike_in_metrics {
    container { python_image }

    input:
    tuple val(run_basename), path(xic_json), path(identifications)
    path spike_ins_table

    output:
    path "${xic_json.baseName}-spikeins.csv"

    """
    extract_spike_metrics.py -itrfp_json ${xic_json} -iidentifications ${identifications} -ispikeins ${spike_ins_table} -ocsv ${xic_json.baseName}-spikeins.csv
    """
}
