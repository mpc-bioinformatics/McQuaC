#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.pia_idents = "$PWD/identifications"  // Peptide Identifications by search machine, which is working with PIA
params.pia_analysis_file = "${baseDir}/example_configurations/pia-analysis.json" // Parameters file to configure PIA (e.g. FDR-Calculation)

// Optional Parameters TODO
params.pia_memory = "8g"  // memory in Java format, passed like "-Xmx[PIA_MEMORY]"
params.pia_outdir = "$PWD/results"  // Output-Directory of the PIA results.


workflow {
    // Run PIA protein inference
    idents = Channel.fromPath(params.pia_idents + "/*.idXML")

    execute_pia(idents)
}


workflow execute_pia {
    take:
        idents
    main:
        // Run PIA protein inference
        pia_compilation(idents)
        pia_analysis(pia_compilation.out)
        pia_extraction(pia_analysis.out)
    emit: 
        pia_analysis.out
        pia_extraction.out
}


process pia_compilation {
    container 'quay.io/biocontainers/pia:1.4.8--hdfd78af_0'
    
    publishDir "${params.pia_outdir}/", mode: "copy"

    input:
    path idXML

    output:
    tuple path("${idXML.baseName}_____pia-compilation.xml"), val("${idXML.baseName}")

    """
    echo "Starting Compilation"
    echo "${idXML}"
    pia -Xmx${params.pia_memory} --compile -o "${idXML.baseName}_____pia-compilation.xml" "${idXML}"
    """
}


process pia_analysis {
    container 'quay.io/biocontainers/pia:1.4.8--hdfd78af_0'

    // Running an analysis with a parameter file
    // The command line allows you to execute an analysis via prior defined analysis in JSON format. 
    // Additionally to the json file, the prior compiled intermediate file must be given.

    publishDir "${params.pia_outdir}/", mode: "copy"

    input:
    tuple path(compilation), val(basename)

    output:
    tuple path("${compilation.simpleName}-piaExport-PSM.mzTab"), path("${compilation.simpleName}-piaExport-peptides.csv"), path("${compilation.simpleName}-piaExport-proteins.mzTab"), val(basename)

    script:
    """
    # TODO: adjust decoy regex in the params file!!
    cat ${params.pia_analysis_file} \
        | sed -e 's;"psmExportFile": "/tmp/piaExport-PSMs.mzTab";"psmExportFile": "${compilation.simpleName}-piaExport-PSM.mzTab";g' \
        | sed -e 's;"peptideExportFile": "/tmp/piaExport-peptides.csv";"peptideExportFile": "${compilation.simpleName}-piaExport-peptides.csv";g' \
        | sed -e 's;"proteinExportFile": "/tmp/piaExport-proteins.mzTab";"proteinExportFile": "${compilation.simpleName}-piaExport-proteins.mzTab";g' > parameters.json
 
    pia -Xmx${params.pia_memory} parameters.json ${compilation}

    # Generate Output, in case of empty results.
    touch ${compilation.simpleName}-piaExport-PSM.mzTab
    touch ${compilation.simpleName}-piaExport-peptides.csv
    touch ${compilation.simpleName}-piaExport-proteins.mzTab
    """
}


process pia_extraction {
    container 'mpc/nextqcflow-python:latest'

    publishDir "${params.pia_outdir}/", mode: "copy"

    input:
    tuple path(psms), path(peptides), path(proteins), val(basename)

    output:
    path "${basename}_____pia_extraction.csv"

    script:
    """
    extract_from_pia_output.py --pia_PSMs $psms --pia_peptides $peptides --pia_proteins $proteins --output ${basename}_____pia_extraction.csv
    """
}
