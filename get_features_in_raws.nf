#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.thermo_raws = "$PWD/raws"  // Folder of Thermo-RAW-files
params.ident_files = "$PWD/raws"  // Folder of mzTab Identification files (already FDR-filtered). Names should be identical to raw files for matching

// Optional Parameters
// Parameters for Feature Detection
// params.resulolution_featurefinder = "-algortihm:mass_trace:mz_tolerance 0.02 -algorithm:isotopic_pattern:mz_tolerance 0.04"  // Parameters for Low Resolution Machines (Q-TOF)
params.resolution_featurefinder = "-algorithm:mass_trace:mz_tolerance 0.004 -algorithm:isotopic_pattern:mz_tolerance 0.005"   // Parameters for High Resolution Machines (LTQ-OrbiTrap)
params.considered_charges_low = "1"  // Charges for the feature finder to use to extract features. In QC this was set to 2:5
params.considered_charges_high = "5"  // Charges for the feature finder to use to extract features. In QC this was set to 2:5

// Output Directory
params.outdir = "$PWD/results"  // Output-Directory of the mzMLs. Here it is <Input_file>.mzML

params.num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)


workflow {
    rawfiles = Channel.fromPath(params.spk_thermo_raws + "/*.raw")
    mztabfiles = Channel.fromPath(params.spk_identification_files + "/*.mzTab")
    retrieve_spikeins(rawfiles, mztabfiles)
}

workflow get_features {
    take:
        rawfiles
        mztabfiles
    main:
        mztabs_tuple = mztabfiles.map {
            file -> tuple(file, file.baseName.split("_____")[0])
        }
        // Match files according to their baseName
        create_baseName_for_raws(rawfiles)
        raw_and_id = create_baseName_for_raws.out.join(
            mztabs_tuple,
            by: 1
        )
        // Convert the file to mzML, where MS1 (peak picked) (for feature finding)
        convert_raw_via_thermorawfileparser(raw_and_id)

        // Get features with OpenMS' feature finder
        run_feature_finder(convert_raw_via_thermorawfileparser.out)

        // Map Identification with Features (using mzTab and featureXML)
        map_features_with_idents(run_feature_finder.out)

        // Retrieve the actual data and report a csv file
        get_statistics_from_featurexml(map_features_with_idents.out)
    emit:
        get_statistics_from_featurexml.out
}

// Stubs for an easy conversion of a channel into a tuple
process create_baseName_for_raws {
    input:
    file raw

    output:
    tuple file(raw), val("$raw.baseName")

    """
    """
}

process convert_raw_via_thermorawfileparser {
    maxForks params.num_procs_conversion
    stageInMode "copy"

    input:
    tuple val(file_base_name), file(raw), file(ident)

    output:
    tuple file("${raw.baseName}.mzML"), file(ident)

    """
    run_thermorawfileparser.sh --format=1 --output_file=${raw.baseName}.mzML --input=${raw} 
    rm ${raw}
    """
}

process run_feature_finder {
    maxForks 3  // TODO Please remove limit, my laptop is not that powerful! 
    stageInMode "copy"

    input:
    tuple file(mzml), file(ident)

    output:
    tuple file("${mzml.baseName}.featureXML"), file(ident)

    """
    run_featurefindercentroided.sh -in ${mzml} -out ${mzml.baseName}.featureXML -algorithm:isotopic_pattern:charge_low ${params.considered_charges_low} -algorithm:isotopic_pattern:charge_high ${params.considered_charges_high} ${params.resolution_featurefinder}
    # run_featurefindermultiplex.sh -in ${mzml} -out ${mzml.baseName}.featureXML
    # We do not use multiplex, it seems to be broken, mem usage way over 40 GB per RAW file failing by "std::bad_alloc"
    rm ${mzml}
    """
}


process map_features_with_idents {
    stageInMode "copy"

    input:
    tuple file(featurexml), file(ident)

    output:
    tuple file("${featurexml.baseName}_with_idents.featureXML"), val("${featurexml.baseName}")
    """
    convert_mztab_to_idxml.py -mztab ${ident} -out_idxml ${ident.baseName}.idXML
    run_idmapper.sh -id ${ident.baseName}.idXML -in ${featurexml} -out ${featurexml.baseName}_with_idents.featureXML
    rm ${featurexml}
    """
}


process get_statistics_from_featurexml {
    stageInMode "copy"

    publishDir "${params.outdir}/", mode:'copy'

    input:
    tuple file(featurexml), val(file_base_name)

    output:
    tuple file(featurexml), file("${file_base_name}_features.csv")

    """
    extract_from_featurexml.py -featurexml ${featurexml} -out_csv ${file_base_name}_features.csv -report_up_to_charge ${params.considered_charges_high}
    """
}
