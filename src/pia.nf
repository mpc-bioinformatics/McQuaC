#!/usr/bin/env nextflow


/**
 * Workflows and processes for protein inference using PIA
 */


nextflow.enable.dsl=2

pia_image = 'quay.io/biocontainers/pia:1.5.5--hdfd78af_0'

params.pia_gb_ram = 16
params.pia_threads = 8

pia_mem_params = "-Xms2g -Xmx" + params.pia_gb_ram + "g"  // additional PIA parameters to set the memory consumption

// Note for CPUs: PIA needs/takes very shortly all available processors, but rather idles on them later.
// On the memory-side, it uses a lot the more PSMs are found...
// Use this information to adjust number of parallel threads.

/*
 * Performs PIA's FDR and protein iference on the given files and levels.
 * 
 * @param identifications search engine results
 * @param do_psm_export whether PSM level data should be returned
 * @param do_peptide_export whether peptide level data should be returned
 * @param do_protein_export whether protein level data should be returned
 * @param fdr_filter whether FDR filtering should be performed (on any level)
 *
 * @return return_files tuples containing the PSM, peptide and protein level results each (may be empty, if level was not returned)
 */
workflow pia_analysis {
    take:
        identifications
        do_psm_export
        do_peptide_export
        do_protein_export
        fdr_filter
    
    main:
        pia_xmls = compile_pia_xmls(identifications)
        analysis_json = prepare_analysis_json(do_psm_export, do_peptide_export, do_protein_export, fdr_filter)
        pia_all_report_files = pia_run_analysis(pia_xmls, analysis_json)

        return_files = pia_all_report_files.psms.collect()
            .concat(pia_all_report_files.peptides.collect())
            .concat(pia_all_report_files.proteins.collect())
            .toList()
            .transpose()
    
    emit:
        return_files
}


/*
 * Performs PIA's FDR and protein iference on the given files and all levels.
 *
 * @param identifications search engine results
 *
 * @return return_files tuples containing the PSM, peptide and protein level results each (may be empty, if level was not returned)
 */
workflow pia_analysis_full {
    take:
        identifications

    main:
        pia_report_files = pia_analysis(identifications, true, true, true, true)
    
    emit:
        pia_report_files
}

/*
 * Performs PIA's FDR and protein iference on the given files and PSM level only.
 *
 * @param identifications search engine results
 * @param fdr_filter whether FDR-filtering should be performed
 *
 * @return psm_file the results on PSM level
 */
workflow pia_analysis_psm_only {
    take:
        identifications
        fdr_filter
    
    main:
        pia_report_files = pia_analysis(identifications, true, false, false, fdr_filter)

    emit:
        pia_report_files
            .toList()
            .transpose()
            .first()
            .flatten()
}

/**
 * Extracts the QC metrics from the PIA results and writes them into a HFD5
 *
 * @param pia_results the PIA results as tuples, as returned by the pia_analysis
 *
 * @return extracted hdf5 metrics written into the HFD5 file
 */
workflow pia_extract_metrics {
    take:
        pia_results

    main:
        extract_csv = pia_extract_csv(pia_results)
    
    emit:
        extract_csv
}


/**
 * Compiles the input files into PIA intermediate files
 *
 * @param identifications search engine results in a PIA usable format
 *
 * @return pia.xml the compilation as pia.xml
 */
process compile_pia_xmls {
    container { pia_image }
    
    cpus {params.pia_threads}
    memory {params.pia_gb_ram + " GB"}

    input:
    path identifications

    output:
    path "${identifications.baseName}.pia.xml"

    """
    pia ${pia_mem_params} --threads ${params.pia_threads} --compile -o ${identifications.baseName}.pia.xml ${identifications}
    """
}

/**
 * Creates a PIA json analysis file with given export parameters
 *
 * @param psm_export true/false either export to PSM level or not
 * @param peptide_export true/false either export to peptide level or not
 * @param protein_export true/false either export to protein level or not
 *
 * @return an pia_analysis.json with defined parameters for the QC
 */
process prepare_analysis_json {
    container { pia_image }

    input:
    val psm_export
    val peptide_export
    val protein_export
    val fdr_filter

    output:
    path "pia_analysis.json"

    """
    pia --example > pia_analysis.json
    
    sed -i 1d pia_analysis.json

    sed -i 's;"createPSMsets": .*,;"createPSMsets": false,;g' pia_analysis.json
    sed -i 's;"psmLevelFileID": .*,;"psmLevelFileID": 1,;g' pia_analysis.json
    sed -i 's;"calculateCombinedFDRScore": .*,;"calculateCombinedFDRScore": false,;g' pia_analysis.json
    if [ ${psm_export} = true ];
    then
      sed -i 's;"psmExportFile":.*,;"psmExportFile": "piaExport-PSMs.mzTab",;g' pia_analysis.json
    else
      sed -i 's;"psmExportFile":.*;;g' pia_analysis.json
    fi

    if [ ${peptide_export} = true ];
    then
      sed -i 's;"inferePeptides":.*,;"inferePeptides": true,;g' pia_analysis.json
      sed -i 's;"peptideExportFile":.*,;"peptideExportFile": "piaExport-peptides.csv",;g' pia_analysis.json
    else
      sed -i 's;"inferePeptides":.*,;"inferePeptides": false,;g' pia_analysis.json
      sed -i 's;"peptideExportFile":.*,;;g' pia_analysis.json
    fi
    sed -i 's;"peptideLevelFileID":.*,;"peptideLevelFileID": 1,;g' pia_analysis.json

    if [ ${protein_export} = true ];
    then
      sed -i 's;"infereProteins":.*,;"infereProteins": true,;g' pia_analysis.json
      sed -i 's;"inferenceMethod":.*,;"inferenceMethod": "inference_spectrum_extractor",;g' pia_analysis.json
      sed -i '/inferenceFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= 0.01"/g;}' pia_analysis.json
      sed -i 's;"scoringBaseScore":.*,;"scoringBaseScore": "psm_fdr_score",;g' pia_analysis.json
      sed -i 's;"scoringPSMs":.*,;"scoringPSMs": "best",;g' pia_analysis.json
      sed -i 's;"proteinExportFile":.*,;"proteinExportFile": "piaExport-proteins.mzTab",;g' pia_analysis.json
      sed -i 's;"proteinExportWithPSMs":.*,;"proteinExportWithPSMs": true,;g' pia_analysis.json
    else
      sed -i 's;"infereProteins":.*,;"infereProteins": false,;g' pia_analysis.json
      sed -i 's;"proteinExportFile":.*,;;g' pia_analysis.json
    fi

    if [ ${fdr_filter} = true ];
    then
      sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= 0.01",\\n    "psm_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
      sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= 0.01",\\n    "peptide_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
      sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= 0.01",\\n    "protein_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
    else
      sed -i '/psmFilters/{n;s/.*//g;}' pia_analysis.json
      sed -i '/peptideFilters/{n;s/.*//g;}' pia_analysis.json
      sed -i '/proteinFilters/{n;s/.*//g;}' pia_analysis.json
    fi
    """
}


/*
 * Runs a PIA analysis with a parameter file
 * The command line allows you to execute an analysis via prior defined analysis in JSON format. 
 * Additionally to the json file, the prior compiled intermediate file must be given.
 *
 * @param pia_xml the pre-compiled PIA XML file
 * @param pia_analysis_file the pre-geenrated PIA json analysis file
 *
 * @return exports for the PSM, peptide and protein level (might be emmpty files, if teh analysis is given like it)
 *
 */
process pia_run_analysis {
    container { pia_image }

    cpus {params.pia_threads}
    memory {params.pia_gb_ram + " GB"}

    input:
    path pia_xml
    path pia_analysis_file

    output:
    path "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-PSM.mzTab", emit: psms
    path "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-peptides.csv", emit: peptides
    path "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-proteins.mzTab", emit: proteins

    """
    filebase="${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}"
    jsonfile="\${filebase}-analysis.json"

    cp ${pia_analysis_file} \${jsonfile}

    touch \${filebase}-piaExport-PSM.mzTab
    touch \${filebase}-piaExport-peptides.csv
    touch \${filebase}-piaExport-proteins.mzTab

    sed -i 's;"psmExportFile": .*;"psmExportFile": "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-PSM.mzTab",;g' \${jsonfile}
    sed -i 's;"peptideExportFile": .*;"peptideExportFile": "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-peptides.csv",;g' \${jsonfile}
    sed -i 's;"proteinExportFile": .*;"proteinExportFile": "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-proteins.mzTab",;g' \${jsonfile}

    pia ${pia_mem_params} --threads ${params.pia_threads} \${jsonfile} ${pia_xml}
    """
}


process pia_extract_csv {
    label 'mcquac_image'

    input:
    tuple path(psm_results), path(peptide_results), path(protein_results)

    output:
    path "${psm_results.name.take(psm_results.name.lastIndexOf('-piaExport-PSM.mzTab'))}-pia_extraction.hdf5"

    script:
    """
    outfile="${psm_results.name.take(psm_results.name.lastIndexOf('-piaExport-PSM.mzTab'))}-pia_extraction.hdf5"
    extract_from_pia_output.py --pia_PSMs ${psm_results} --pia_peptides ${peptide_results} --pia_proteins ${protein_results} --out_hdf5 \${outfile}
    """
}
