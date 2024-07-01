#!/usr/bin/env nextflow
nextflow.enable.dsl=2

python_image = 'mpc/nextqcflow-python:latest'

params.openms_threads = 8   // hardcoded for now, number of threads used by OpenMS
params.min_charge = 2       // hardcoded for now
params.max_charge = 5       // hardcoded for now

/**
 * Extracts features from peak picked MS1 spectra
 * 
 * @param mzmls Channel of mzML files
 * @param mztabfiles Channel of mzTab files
 * @param feature_finder_params Channel of feature finder parameters
 */
workflow get_feature_metrics {
    take:
        mzmls  // MS1 should be peak picked for feature-finding
        mztabfiles
        comet_params
    
    main:
        // filter out empty spectra and chromatograms
        filtered_mzml = filter_mzml(mzmls)

        // Get features with OpenMS feature finder
        feature_finder_params = get_feature_finder_params_from_comet_params(comet_params)
        feature_xml = run_feature_finder(filtered_mzml, feature_finder_params)

        // map identification to features (using mzTab and featureXML)
        runbase_to_featurexml_and_mztab = 
            feature_xml.map {
                file -> tuple(file.name.take(file.name.lastIndexOf('.filtered.featureXML')), file)
            }.join(
                mztabfiles.map{
                    file -> tuple(file.name.take(file.name.lastIndexOf('-piaExport-PSM.mzTab')), file)
                },
                by: 0
            )
        feature_xml_identified = map_features_to_idents(runbase_to_featurexml_and_mztab)

        // Retrieve the actual data and report a csv file
        feature_metrics = get_metrics_from_featurexml(feature_xml_identified)
        
    emit:
         feature_metrics
}

/**
 * This process is used to convert the comet parameters into the parameters needed for the feature finder.
 * 
 * @param comet_params The comet parameters file
 * @return The feature finder parameters (value channel)
 */
process get_feature_finder_params_from_comet_params {
	container { python_image }

	input:
	path comet_params

	output:
	stdout

	"""
	comet_params_to_feature_finder_params.py -c ${comet_params}
	"""
}

process run_feature_finder {
    container { python_image }

    cpus { params.openms_threads }

    input:
    path mzml
    val feature_finder_params

    output:
    path "${mzml.baseName}.featureXML"

    """
    FeatureFinderMultiplex -in ${mzml} -out ${mzml.baseName}.featureXML -threads ${params.openms_threads} \
        -algorithm:labels "" \
        -algorithm:charge "${params.min_charge}:${params.max_charge}" \
        -algorithm:spectrum_type centroid \
        ${feature_finder_params}
    """
}

process map_features_to_idents {
    container { python_image }

    cpus { params.openms_threads }

    input:
    tuple val(runBaseName), path(featurexml), path(ident)

    output:
    path "${featurexml.baseName}_with_idents.featureXML"

    """
    convert_mztab_to_idxml.py -mztab ${ident} -out_idxml ${ident.baseName}.idXML
    IDMapper -in ${featurexml} -out ${featurexml.baseName}_with_idents.featureXML -threads ${params.openms_threads} \
        -id ${ident.baseName}.idXML 
    """
}

process get_metrics_from_featurexml {
    container { python_image }

    input:
    path featurexml

    output:
    path "${featurexml.name.take(featurexml.name.lastIndexOf('.filtered_with_idents.featureXML'))}-features.csv"

    """
    extract_from_featurexml.py -featurexml ${featurexml} \
        -report_up_to_charge ${params.max_charge} \
        -out_csv ${featurexml.name.take(featurexml.name.lastIndexOf('.filtered_with_idents.featureXML'))}-features.csv 
    """
}

/**
 * Removes chromatrograms and empty peak lists from mzML
 * Necessary to prevent memory issues with FeatureFinder
 */
process filter_mzml {
    container { python_image }

    cpus { params.openms_threads }

    input:
    path(mzml)

    output:
    path "${mzml.baseName}.filtered.mzML"

    """
    # Filter mzML for MS1 spectra
    FileFilter -in ${mzml} -out ${mzml.baseName}.filtered.mzML -threads ${params.openms_threads} \
        -peak_options:remove_chromatograms \
        -peak_options:remove_empty \
        -peak_options:sort_peaks \
        -peak_options:zlib_compression true
    """
}