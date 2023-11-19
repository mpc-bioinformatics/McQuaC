#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.pia_idents = "$PWD/identifications"  // Peptide Identifications by search machine, which is working with PIA

// Optional Parameters TODO
params.pia_outdir = "$PWD/results"  // Output-Directory of the PIA results.
params.pia_parameters_file = ""  // Parameters file to configure PIA (e.g. FDR-Calculation)

params.pia_memory = "8g"  // memory in Java format, passed like "-Xmx[PIA_MEMORY]"
params.pia_executable = "${baseDir}/bin/pia/pia-1.4.8.jar"  // Executable to the pia jar
params.pia_analysis_file = "${baseDir}/example_configurations/pia-analysis.json"

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
    java -Xmx${params.pia_memory} -jar "${params.pia_executable}" --compile -o "${idXML.baseName}_____pia-compilation.xml" "${idXML}"
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
    # TODO: adjust decoy regex in the params file!!
    cat ${params.pia_analysis_file} \
        | sed -e 's;"psmExportFile": "/tmp/piaExport-PSMs.mzTab";"psmExportFile": "${compilation.simpleName}-piaExport-PSM.mzTab";g' \
        | sed -e 's;"peptideExportFile": "/tmp/piaExport-peptides.csv";"peptideExportFile": "${compilation.simpleName}-piaExport-peptides.csv";g' \
        | sed -e 's;"proteinExportFile": "/tmp/piaExport-proteins.mzTab";"proteinExportFile": "${compilation.simpleName}-piaExport-proteins.mzTab";g' > parameters.json
 
    java -Xmx${params.pia_memory} -jar "${params.pia_executable}" parameters.json ${compilation}

    # Generate Output, in case of empty results.
    touch ${compilation.simpleName}-piaExport-PSM.mzTab
    touch ${compilation.simpleName}-piaExport-peptides.csv
    touch ${compilation.simpleName}-piaExport-proteins.mzTab
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
}
