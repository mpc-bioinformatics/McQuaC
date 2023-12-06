#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.spk_raw_spectra = "$PWD/raws"  // Databasefile in SP-EMBL
params.spk_identification_files = "$PWD/idents" // Identification-Files in mzTAB-Format already FDR-filtered (e.g. by PIA)

// Optional Parameters
params.spk_spike_ins = "${baseDir}/example_configurations/spike_ins.csv" // The Spike-Ins we also added in the Identification. Defaults to MPC-SPIKEINS
params.spk_outdir = "$PWD/results"  // Output-Directory of the XICs (and corrected retention-time XICs with the help of the identification) as a csv-file. Here it is <Input-Raw-File>_spikeins.csv
params.spk_num_procs_extraction = Runtime.runtime.availableProcessors()  // Number of process used to extract (CAUTION: This can be very resource intensive!)


workflow {
    rawfiles = Channel.fromPath(params.spk_raw_spectra + "/*.{raw,d}", type: "any")
    ident_files = Channel.fromPath(params.spk_identification_files + "/*.mzTab")
    retrieve_spikeins(rawfiles, ident_files)
}

workflow retrieve_spikeins {
    take:
        raw_files
        mztabfiles
    main:
        // Create a value channel for concurrent processing
        // TODO: should be a parameter of the workflow to support other spike-ins when calling from main
        spikeins = Channel.fromPath(params.spk_spike_ins).first()

        // Map the raw files to their corresponding mzTab files
        // Format: (basename, raw, mztab)
        rawfiles_tuple = raw_files.map {
            file -> tuple(file.baseName, file)
        }
        mztabs_tuple = mztabfiles.map {
            file -> tuple(file.baseName.split("_____")[0], file)
        }
        raw_and_id = rawfiles_tuple.join(
            mztabs_tuple,
            by: 0
        )

        // Finally, we generate the input json, retrieve it via trfp and parse back this results into a csv-format
        json_and_association = generate_json_and_association(raw_and_id, spikeins)


        json_and_association.branch {
            thermo: it[0].getExtension() == 'raw'
            bruker: it[0].getExtension() == 'd'
        }.set{ filtered_json_and_association }


        thermo_xics = retrieve_xics_from_thermo_raw_spectra(filtered_json_and_association.thermo)
        bruker_xics = retrieve_xics_from_bruker_raw_spectra(filtered_json_and_association.bruker)

        xics = thermo_xics.concat(bruker_xics)

        get_statistics(xics)

    emit:
        get_statistics.out
}

// Uses the Identification file and Spike Ins and maps according to the accession identified spike ins. We also generate here the query for TRFP 
process generate_json_and_association {
    container 'mpc/nextqcflow-python:latest'

    input:
    tuple val(file_base_name), path(raw), path(ident)
    path spike_ins

    output:
    tuple path(raw), path("trfp_input.json"), path("association.txt")

    """
    xics_to_json.py -icsv $spike_ins -iidents $ident -ojson trfp_input.json -oassociation association.txt
    """
}

// Actual retrieval of the XICs using TRFP (Wrapper)
process retrieve_xics_from_thermo_raw_spectra {
    container 'quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0'

    maxForks params.spk_num_procs_extraction

    input:
    tuple path(raw), path(trfp_input), path(associations)

    output:
    tuple path(associations), path("${raw.baseName}.json")

    """
    thermorawfileparser xic -i $raw -j $trfp_input 
    """
}

// Actual retrieval of the XICs using TRFP (Wrapper)
process retrieve_xics_from_bruker_raw_spectra {
    container 'mpc/nextqcflow-python:latest'

    maxForks params.spk_num_procs_extraction

    input:
    tuple path(raw), path(trfp_input), path(associations)

    output:
    tuple path(associations), path("${raw.baseName}.json")

    """
    extract_xic_bruker.py -d_folder ${raw} -in_json ${trfp_input} -out_json ${raw.baseName}.json
    """
}

// Parsing back the results from TRFP (in JSON) into a CSV-format (while also considering identifications)
process get_statistics {
    publishDir "${params.spk_outdir}/", mode:'copy'

    input:
    tuple val(associations), path(trfp_spike_ins_json)

    output:
    path("${trfp_spike_ins_json.baseName}_____spikeins.csv")

    """
    trfp_json_to_table.py -itrfp_json $trfp_spike_ins_json -iassociations $associations -ocsv ${trfp_spike_ins_json.baseName}_____spikeins.csv
    """
}
