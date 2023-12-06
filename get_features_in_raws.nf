#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.gf_thermo_raws = "$PWD/raws"  // Folder of Thermo-RAW-files
params.gf_ident_files = "$PWD/raws"  // Folder of mzTab Identification files (already FDR-filtered). Names should be identical to raw files for matching
params.feature_finder_params = "-algorithm:mz_tolerance 5.0 -algorithm:mz_unit ppm" // Mass spectrometer specific parameters for FeatureFinder

// Optional Parameters
// Parameters for Feature Detection
// params.resulolution_featurefinder = "-algortihm:mass_trace:mz_tolerance 0.02 -algorithm:isotopic_pattern:mz_tolerance 0.04"  // Parameters for Low Resolution Machines (Q-TOF)
params.gf_resolution_featurefinder = "-algorithm:mass_trace:mz_tolerance 0.004 -algorithm:isotopic_pattern:mz_tolerance 0.005"   // Parameters for High Resolution Machines (LTQ-OrbiTrap)
params.gf_considered_charges_low = "1"  // Charges for the feature finder to use to extract features. In QC this was set to 2:5
params.gf_considered_charges_high = "6"  // Charges for the feature finder to use to extract features. In QC this was set to 2:5
params.additional_dinosaur_settings = "" // Additional Parameters which can be set for Dinosaur


// Output Directory
params.gf_outdir = "$PWD/results"  // Output-Directory of the mzMLs. Here it is <Input_file>.mzML

params.gf_num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)


workflow {
    mzmls = Channel.fromPath(params.gf_thermo_raws + "/*.mzML")
    mztabfiles = Channel.fromPath(params.gf_ident_files + "/*.mzTab")
    feature_finder_params = Channel.value(params.feature_finder_params)
    get_features(mzmls, mztabfiles, feature_finder_params)
}

/**
 * Extracts features from peak picked MS1 spectra
 * 
 * @param mzmls Channel of mzML files
 * @param mztabfiles Channel of mzTab files
 * @param feature_finder_params Channel of feature finder parameters
 */
workflow get_features {
    take:
        mzmls  // MS1 should be peak picked for feature-finding
        mztabfiles
        feature_finder_params
    main:
        // Match files according to their baseName
        mztabs_tuple = mztabfiles.map {
            file -> tuple(file, file.baseName.split("_____")[0])
        }
        mzmls_tuple = mzmls.map {
            file -> tuple(file, file.baseName)
        }
        mzml_and_id = mzmls_tuple.join(
            mztabs_tuple,
            by: 1
        ).map {
            it -> tuple(it[1], it[2])
        }

        // Get features with OpenMS' or Dinosaur feature finder
        run_feature_finder(mzml_and_id, feature_finder_params)

        // Map Identification with Features (using mzTab and featureXML)
        map_features_with_idents(run_feature_finder.out)

        // Retrieve the actual data and report a csv file
        get_statistics_from_featurexml(map_features_with_idents.out)
    emit:
        get_statistics_from_featurexml.out
}

process run_feature_finder {
    container 'mpc/nextqcflow-python:latest'

    input:
    tuple path(mzml), path(ident)
    val feature_finder_params

    output:
    tuple path("${mzml.baseName}.featureXML"), path("${mzml.baseName}.hills.csv"), path(ident)

    """
    # Centroided FF
    # FeatureFinderCentroided -in ${mzml} -out ${mzml.baseName}.featureXML -algorithm:isotopic_pattern:charge_low ${params.gf_considered_charges_low} -algorithm:isotopic_pattern:charge_high ${params.gf_considered_charges_high} ${params.gf_resolution_featurefinder}
    # touch ${mzml.baseName}.hills.csv

    # Multiplex FF
    # Suggested by OpenMS developers, even if there is no multiplexing (just pass empty labels)
    # Ensure params fit to your mass spectrometer, otherise this will eat up your memory faster than Chrome.
    # Just use charge state 2 - 5 as this are common charge states for peptides lower won't be identified anyway
    FeatureFinderMultiplex -in ${mzml} -out ${mzml.baseName}.featureXML \
        -algorithm:labels "" \
        -algorithm:charge "2:5" \
        -threads ${params.gf_num_procs_conversion} \
        ${params.feature_finder_params}
    touch ${mzml.baseName}.hills.csv
    
    # Dinosaur FF
    # java -jar \$(get_cur_bin_dir.sh)/Dinosaur.jar \
    #    --writeHills \
    #    --writeMsInspect \
    #    ${params.additional_dinosaur_settings} --mzML ${mzml}
    # FileConverter -in ${mzml.baseName}.msInspect.tsv -out ${mzml.baseName}.featureXML
    """
}

process map_features_with_idents {
    container 'mpc/nextqcflow-python:latest'

    input:
    tuple path(featurexml), path(hills), path(ident)

    output:
    tuple path("${featurexml.baseName}_with_idents.featureXML"), path(hills), val("${featurexml.baseName}")
    """
    convert_mztab_to_idxml.py -mztab ${ident} -out_idxml ${ident.baseName}.idXML
    IDMapper -id ${ident.baseName}.idXML -in ${featurexml} -out ${featurexml.baseName}_with_idents.featureXML
    """
}

process get_statistics_from_featurexml {
    container 'mpc/nextqcflow-python:latest'

    publishDir "${params.gf_outdir}/", mode:'copy'

    input:
    tuple path(featurexml), path(hills), val(file_base_name)

    output:
    tuple path(featurexml), path("${file_base_name}_____features.csv")

    """
    extract_from_featurexml.py -featurexml ${featurexml} -hills ${hills} -out_csv ${file_base_name}_____features.csv -report_up_to_charge ${params.gf_considered_charges_high}
    """
}
