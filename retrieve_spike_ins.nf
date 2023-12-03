#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.spk_raw_spectra = "$PWD/raws"  // Databasefile in SP-EMBL
params.spk_identification_files = "$PWD/idents" // Identification-Files in mzTAB-Format already FDR-filtered (e.g. by PIA)

// Optional Parameters
params.spk_spike_ins = "$PWD/example_configurations/spike_ins.csv" // The Spike-Ins we also added in the Identification. Defaults to MPC-SPIKEINS
params.spk_outdir = "$PWD/results"  // Output-Directory of the XICs (and corrected retention-time XICs with the help of the identification) as a csv-file. Here it is <Input-Raw-File>_spikeins.csv
params.spk_num_procs_extraction = Runtime.runtime.availableProcessors()  // Number of process used to extract (CAUTION: This can be very resource intensive!)


workflow {
    rawfiles = Channel.fromPath(params.spk_raw_spectra + "/*.{raw,d}", type: "any")
    ident_files = Channel.fromPath(params.spk_identification_files + "/*.mzTab")
    retrieve_spikeins(rawfiles, ident_files)
}

workflow retrieve_spikeins {
    take:
        rawfiles
        mztabfiles
    main:
        spikeins = Channel.from(file(params.spk_spike_ins))
        
        // Match files according to their baseName
        rawfiles_tuple = rawfiles.map {
            file -> tuple(file.baseName, file)
        }
        mztabs_tuple = mztabfiles.map {
            file -> tuple(file.baseName.split("_____")[0], file)
        }
        raw_and_id = rawfiles_tuple.join(
            mztabs_tuple,
            by: 0
        )

        // And combine all with the spike ins 
        raw_id_spikes = raw_and_id.combine(spikeins)

        // Finally, we generate the input json, retrieve it via trfp and parse back this results into a csv-format
        generate_json_and_association(raw_id_spikes)
        retrieve_xics_from_raw_spectra(generate_json_and_association.out)
        get_statistics(retrieve_xics_from_raw_spectra.out)

    emit:
        get_statistics.out
}

// Uses the Identification file and Spike Ins and maps according to the accession identified spike ins. We also generate here the query for TRFP 
process generate_json_and_association {
    input:
    tuple val(file_base_name), file(raw), file(ident), file(spike_ins)

    output:
    tuple file(raw), file("trfp_input.json"), file("association.txt")

    """
    xics_to_json.py -icsv $spike_ins -iidents $ident -ojson trfp_input.json -oassociation association.txt
    """
}

// Actual retrieval of the XICs using TRFP (Wrapper)
process retrieve_xics_from_raw_spectra {
    maxForks params.spk_num_procs_extraction
    stageInMode "copy"

    input:
    tuple file(raw), file(trfp_input), file(associations)

    output:
    tuple file(associations), file("${raw.baseName}.json")


    """
    if [[ "${raw}" == *.raw ]]; then
        thermorawfileparser xic -i $raw -j $trfp_input 
    fi

    if [[ "${raw}" == *.d ]]; then
        extract_xic_bruker.py -d_folder ${raw} -in_json ${trfp_input} -out_json ${raw.baseName}.json
    fi 
    """
}

// Parsing back the results from TRFP (in JSON) into a CSV-format (while also considering identifications)
process get_statistics {
    publishDir "${params.spk_outdir}/", mode:'copy'

    input:
    tuple val(associations), file(trfp_spike_ins_json)

    output:
    file("${trfp_spike_ins_json.baseName}_____spikeins.csv")

    """
    trfp_json_to_table.py -itrfp_json $trfp_spike_ins_json -iassociations $associations -ocsv ${trfp_spike_ins_json.baseName}_____spikeins.csv
    """
}
