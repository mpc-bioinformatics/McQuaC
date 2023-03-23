#!/usr/bin/env nextflow
/* here stands super valuable information

Please use camel case for process names and snake case for variable names. Go all caps on global, non-dynamic variables. 
Have fun!

Example run for PIA: nextflow run pia.nf --idents /home/maike/Programmierprojekte/QC-Nextflow/pia-tutorial/data/identifications/ --outdir /home/maike/Programmierprojekte/QC-Nextflow/result/
*/

PROJECT_DIR = workflow.projectDir
NFX_WORK = PROJECT_DIR + "work/"
println(PROJECT_DIR)
params.help = false

if(params.help) {
	println(""" \nusage : ~/nextflow main.nf [operators] --input inputFiles
		example : ~/nextflow main.nf --input \"path/fastas/\"
		operators: --help calls this help option\n """)
	System.exit(0)
	} else {
		dir= workflow.projectDir.getParent() + "/results/"
		

//rawfiles_pia = Channel.fromPath(PROJECT_DIR  + "/results/mzid/*.mzid")

workflow {
	converter()
	//comet(converter.out)
	msgfplus(converter.out)
	//msgfplus.out.view()
    pia(msgfplus.out)
}
//add parameters to workflow : include {sayHello} from './some/module' addParams(foo: 'Ciao')
include {convert_raw_via_thermorawfileparser} from PROJECT_DIR + '/convert_to_mgf_thermorawfileparser.nf'

workflow converter{
    // Convert the file to MGF
	// Needs to have mono installed on linux system (mono-complete)
	// required parameters: params.thermo_raws
    rawfiles = Channel.fromPath(params.thermo_raws + "/*.raw")
	main:
    	convert_raw_via_thermorawfileparser(rawfiles)
	emit:
		convert_raw_via_thermorawfileparser.out
}

params.search_parameter_file = "$PWD/example_configurations/comet_config.txt" //Search Parameters for Comet
include {comet_search_mgf} from PROJECT_DIR + "/identification_via_comet.nf"
workflow comet{
    // Get all MGF files which should be identified
 //   mgfs = Channel.fromPath(converter.out)

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.search_parameter_file)


	take: data
	main: 
    // Start search
		combined_channel = fasta_file
        	.combine(modifications_file)
        	.combine(data)
		combined_channel.view()
    	comet_search_mgf(combined_channel)
	emit:
		comet_search_mgf.out
}

include {msgfplus_buildsa; msgfplus_search_mgf} from PROJECT_DIR + "/identification_via_msgfplus.nf"
workflow msgfplus{
    // Get all MGF files which should be identified
 //   mgfs = Channel.fromPath(converter.out)

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.search_parameter_file)

	take: data

	main: 
    // Build indexed fasta for MSGFPLUS
    	msgfplus_buildsa(fasta_file)
    // Combined channel replicated the indexed fasta for each MGF to be reused
    	combined_channel = fasta_file
        	.combine(modifications_file)
        	.combine(data)
        	.combine(msgfplus_buildsa.out.toList())
	// Start search
    	msgfplus_search_mgf(combined_channel)
		msgfplus_search_mgf.out.view()

	emit:
		msgfplus_search_mgf.out
}


include {pia_compilation; pia_analysis} from PROJECT_DIR + '/pia.nf'

workflow pia{
    // Run PIA protein inference
	//params.idents = "/home/maike/Programmierprojekte/QC-Nextflow/pia-tutorial/data/identifications/"
    
	take: data
	main:
    	pia_compilation(data)
    	pia_analysis(pia_compilation.out)
	emit: 
		pia_analysis.out
}
}
