#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.thermo_raws = "$PWD/raws"  // Databasefile in SP-EMBL
params.identification_files = "$PWD/idents" // Identification-Files in mzTAB-Format already FDR-filtered (e.g. by PIA)
params.spike_ins = "$PWD/parameter_files/spike_ins.csv" // The Spike-Ins we also added in the Identification

// Optional Parameters
params.outdir = "$PWD/results"  // Output-Directory of the XICs (and corrected retention-time XICs with the help of the identification) as a csv-file. Here it is <Input-Raw-File>_spikeins.csv
params.num_procs_extraction = Runtime.runtime.availableProcessors()  // Number of process used to extract (CAUTION: This can be very resource intensive!)

workflow {
    // Get files from Folder
    rawfiles = Channel.fromPath(params.thermo_raws + "/*.raw")
    mztabfiles = Channel.fromPath(params.identification_files + "/*.mzTab")
    
    // Match files according to their baseName
    create_baseName_for_raws(rawfiles)
    create_baseName_for_idents(mztabfiles)
    raw_and_id = create_baseName_for_raws.out.join(
        create_baseName_for_idents.out,
        by: 1
    )

    // And combine all with the spike ins 
    raw_id_spikes = raw_and_id.combine(Channel.from(file(params.spike_ins)))

    // Finally, we generate the input json, retrieve it via trfp and parse back this results into a csv-format
    generate_json_and_association(raw_id_spikes)
    retrieve_via_thermorawfileparser(generate_json_and_association.out)
    get_statistics(retrieve_via_thermorawfileparser.out)
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

process create_baseName_for_idents {
    input:
    file mztab

    output:
    tuple file(mztab), val("$mztab.baseName")

    """
    """
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
process retrieve_via_thermorawfileparser {
    maxForks params.num_procs_extraction
    stageInMode "copy"

    input:
    tuple file(raw), file(trfp_input), file(associations)

    output:
    tuple file(associations), file("${raw.baseName}.json")


    """
    run_thermorawfileparser.sh xic -i $raw -j $trfp_input 
    rm ${raw}
    """
}

// Parsing back the results from TRFP (in JSON) into a CSV-format (while also considering identifications)
process get_statistics {
    publishDir "${params.outdir}/", mode:'copy'

    input:
    tuple val(associations), file(trfp_spike_ins_json)

    output:
    file("${trfp_spike_ins_json.baseName}_spikeins.csv")

    """
    trfp_json_to_table.py -itrfp_json $trfp_spike_ins_json -iassociations $associations -ocsv ${trfp_spike_ins_json.baseName}_spikeins.csv
    """
}
