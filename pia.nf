#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.idents = "$PWD/identifications"  // Peptide Identifications by search machine
params.outdir = "$PWD/results"  // Output-Directory of the PIA results. 3 files are expected.

params.additional_params = ""
params.num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

workflow pia{
    // Run PIA protein inference
    rawfiles = Channel.fromPath(params.idents + "/*.idXML")
    pia_compilation(rawfiles)
    pia_analysis(pia_compilation.out)
}

process pia_compilation {
    // Compiling files into a PIA intermediate file
    maxForks params.num_procs_conversion
    stageInMode "copy"
    
    publishDir "${params.outdir}/", mode:'copy'

    input:
    file idXML

    output:
    file "${idXML.baseName}-pia-compilation.xml"

    """
    echo "Starting Compilation"
    echo "${idXML}"
    java -jar "${baseDir}/bin/pia-1.4.7/pia-1.4.7.jar" --compile -o "${idXML.baseName}-pia-compilation.xml"  ${idXML}
    """
}

process pia_analysis {
    // Running an analysis with a parameter file
    // The command line allows you to execute an analysis via prior defined analysis in JSON format. 
    // Additionally to the json file, the prior compiled intermediate file must be given.
    publishDir "${params.outdir}/", mode:'copy'

    input:
    file compilation

    output:
    file "${compilation.simpleName}-piaExport-PSM.mzTab"
    // file "${compilation.simpleName}-piaExport-peptides.csv"
    // file "${compilation.simpleName}-piaExport--proteins.mzid"


    script:
    """
    cat ${baseDir}/bin/pia-1.4.7/pia-analysis.json \
        | sed -e 's;"psmExportFile": "/tmp/piaExport-PSMs.mzTab";"psmExportFile": "${compilation.simpleName}-piaExport-PSM.mzTab";g' \
        | sed -e 's;"peptideExportFile": "/tmp/piaExport-peptides.csv";"peptideExportFile": "${compilation.simpleName}-piaExport-peptides.csv";g' \
        | sed -e 's;"proteinExportFile": "/tmp/piaExport-proteins.mzid";"proteinExportFile": "${compilation.simpleName}-piaExport--proteins.mzid";g' > parameters.json
 
    java -jar "${baseDir}/bin/pia-1.4.7/pia-1.4.7.jar" parameters.json ${compilation}
    """
}

