#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// this version of the workflow calls PIA as it is distributed by Bioconda.

// Parameters required for standalone execution
params.pia_idents = "$PWD/identifications"  // Peptide Identifications by search machine, which is working with PIA

// Optional Parameters TODO
params.pia_outdir = "$PWD/results"  // Output-Directory of the PIA results.
params.pia_parameters_file = ""  // Parameters file to configure PIA (e.g. FDR-Calculation)

params.pia_memory = "8g"        // memory in Java format, passed like "-Xmx[PIA_MEMORY]"

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
    // Compiling files into a PIA intermediate file
    stageInMode "symlink"
    
    publishDir "${params.pia_outdir}/", mode: "copy"

    input:
    file idXML

    output:
    tuple file("${idXML.baseName}_____pia-compilation.xml"), val("${idXML.baseName}")

    """
    echo "Starting Compilation"
    echo "${idXML}"
    pia -Xmx${params.pia_memory} --compile -o "${idXML.baseName}_____pia-compilation.xml" "${idXML}"
    """
}


process pia_analysis {
    // Running an analysis with a parameter file
    // The command line allows you to execute an analysis via prior defined analysis in JSON format. 
    // Additionally to the json file, the prior compiled intermediate file must be given.
    stageInMode "symlink"

    publishDir "${params.pia_outdir}/", mode: "copy"

    input:
    tuple file(compilation), val(basename)

    output:
    tuple file("${compilation.simpleName}-piaExport-PSM.mzTab"), file("${compilation.simpleName}-piaExport-peptides.csv"), file("${compilation.simpleName}-piaExport-proteins.mzTab"), val(basename)

    script:
    """
    #TODO REMOVE HARD_CODED parameters-file and make it available from outside the script!
    # TODO: adjust decoy regex in the params file!!
    cat ${baseDir}/example_configurations/pia-analysis.json \
        | sed -e 's;"psmExportFile": "/tmp/piaExport-PSMs.mzTab";"psmExportFile": "${compilation.simpleName}-piaExport-PSM.mzTab";g' \
        | sed -e 's;"peptideExportFile": "/tmp/piaExport-peptides.csv";"peptideExportFile": "${compilation.simpleName}-piaExport-peptides.csv";g' \
        | sed -e 's;"proteinExportFile": "/tmp/piaExport-proteins.mzTab";"proteinExportFile": "${compilation.simpleName}-piaExport-proteins.mzTab";g' > parameters.json
 
    pia -Xmx${params.pia_memory} parameters.json ${compilation}
    """
}


process pia_extraction {
    stageInMode "symlink"

    publishDir "${params.pia_outdir}/", mode: "copy"

    input:
    tuple file(psms), file(peptides), file(proteins), val(basename)

    output:
    file "${basename}_____pia_extraction.csv"

    script:
    """
    python ${baseDir}/bin/extract_from_pia_output.py --pia_PSMs $psms --pia_peptides $peptides --pia_proteins $proteins --output ${basename}_____pia_extraction.csv
    """
    //# DEBUGGING
    //#echo "test_col1,test_col2\nvalue1,value2\n" > ${basename}_____pia_extraction.csv
}
