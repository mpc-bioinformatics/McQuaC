#!/usr/bin/env nextflow

// 
// Workflows and processes for protein inference using PIA
// 

// Note for CPUs: PIA needs/takes very shortly all available processors, but rather idles on them later.
// On the memory-side, it uses a lot the more PSMs are found...
// Use this information to adjust number of parallel threads.

/*
 * Performs PIA's FDR and protein iference on the given files and levels.
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
        pia_threads
        pia_gb_ram
        pia_prefilter_threshold

    main:
        pia_xmls = compile_pia_xmls(identifications, pia_threads, pia_gb_ram)

        if (pia_prefilter_threshold > 0) {
            // perform the pre-filtering
            pia_xmls = perform_pia_prefiltering(pia_xmls, pia_prefilter_threshold, pia_threads, pia_gb_ram)
        }

        analysis_json = prepare_analysis_json(do_psm_export, do_peptide_export, do_protein_export, fdr_filter, true, 0.01)
        pia_all_report_files = pia_run_analysis(pia_xmls, analysis_json, pia_threads, pia_gb_ram, "mzTab")

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
 * @return return_files tuples containing the PSM, peptide and protein level results each (may be empty, if level was not returned)
 */
workflow pia_analysis_full {
    take:
        identifications
        pia_threads
        pia_gb_ram
        pia_prefilter_threshold

    main:
        pia_report_files = pia_analysis(
            identifications,
            true,
            true,
            true,
            true,
            pia_threads,
            pia_gb_ram,
            pia_prefilter_threshold
        )
    
    emit:
        pia_report_files
}

/*
 * Performs PIA's FDR and protein iference on the given files and PSM level only.
 *
 * @return psm_file the results on PSM level
 */
workflow pia_analysis_psm_only {
    take:
        identifications
        fdr_filter
        pia_threads
        pia_gb_ram
    
    main:
        pia_report_files = pia_analysis(
            identifications,
            true,
            false,
            false,
            fdr_filter,
            pia_threads,
            pia_gb_ram,
            -1
        )

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


workflow perform_pia_prefiltering {
    take:
        pia_xmls
        pia_prefilter_threshold
        pia_threads
        pia_gb_ram

    main:
        // perform the pre-filtering
        prefilter_analysis_json = prepare_analysis_json(true, false, false, true, false, pia_prefilter_threshold)
        pia_prefiltered = pia_run_analysis(pia_xmls, prefilter_analysis_json, pia_threads, pia_gb_ram, "mzid")
        prefiltered_pia_xmls = compile_pia_xmls(pia_prefiltered.psms, pia_threads, pia_gb_ram)

    emit:
    prefiltered_pia_xmls
}


/**
 * Compiles the input files into PIA intermediate files
 *
 * @return pia.xml the compilation as pia.xml
 */
process compile_pia_xmls {
    label 'pia_image'
    
    cpus {pia_threads}
    memory "${pia_gb_ram}.GB"

    input:
    path identifications
    val pia_threads
    val pia_gb_ram

    output:
    path "${identifications.baseName}.pia.xml"

    script:
    """
    pia -Xms2g -Xmx${pia_gb_ram}g --threads ${pia_threads} --compile -o ${identifications.baseName}.pia.xml ${identifications}
    """
}

/**
 * Creates a PIA json analysis file with given export parameters
 *
 * @return an pia_analysis.json with defined parameters for the QC
 */
process prepare_analysis_json {
    label 'pia_image'

    cpus 1
    memory "4.GB"

    input:
    val psm_export
    val peptide_export
    val protein_export
    val fdr_filter
    val remove_decoys       // remove decoys works only together with fdr_filer == true
    val fdr_threshold

    output:
    path "pia_analysis.json"

    script:
    """
    pia --example > pia_analysis.json
    
    # delete the first row of the file
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
      if [ ${remove_decoys} = true ];
      then
        sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}",\\n    "psm_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
        sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}",\\n    "peptide_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
        sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= ${fdr_threshold}",\\n    "protein_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
      else
        sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
        sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
        sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= ${fdr_threshold}"/g;}' pia_analysis.json
      fi
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
 * @return exports for the PSM, peptide and protein level (might be emmpty files, if teh analysis is given like it)
 *
 */
process pia_run_analysis {
    label 'pia_image'

    cpus {pia_threads}
    memory "${pia_gb_ram}.GB"

    input:
    path pia_xml
    path pia_analysis_file
    val pia_threads
    val pia_gb_ram
    val PSMformat           // needs to be changed between pre-processing (mzid9) and final (mzTab) 

    output:
    path "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-PSM.${PSMformat}", emit: psms
    path "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-peptides.csv", emit: peptides
    path "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-proteins.mzTab", emit: proteins

    script:
    """
    filebase="${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}"
    jsonfile="\${filebase}-analysis.json"

    cp ${pia_analysis_file} \${jsonfile}

    touch \${filebase}-piaExport-PSM.${PSMformat}
    touch \${filebase}-piaExport-peptides.csv
    touch \${filebase}-piaExport-proteins.mzTab

    sed -i 's;"psmExportFile": .*;"psmExportFile": "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-PSM.${PSMformat}",;g' \${jsonfile}
    sed -i 's;"peptideExportFile": .*;"peptideExportFile": "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-peptides.csv",;g' \${jsonfile}
    sed -i 's;"proteinExportFile": .*;"proteinExportFile": "${pia_xml.name.take(pia_xml.name.lastIndexOf('.pia.xml'))}-piaExport-proteins.mzTab",;g' \${jsonfile}

    pia -Xms2g -Xmx${pia_gb_ram}g --threads ${pia_threads} \${jsonfile} ${pia_xml}
    """
}


process pia_extract_csv {
    label 'mcquac_image'

    cpus 1
    memory "8.GB"

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
