#!/usr/bin/env nextflow
/* here stands super valuable information

Please use camel case for process names and snake case for variable names. Go all caps on global, non-dynamic variables. 
Have fun!
*/

PROJECT_DIR = workflow.projectDir
NFX_WORK = PROJECT_DIR + "work/"

params.help = false

if(params.help) {
	println(""" \nusage : ~/nextflow main.nf [operators] --input inputFiles
		example : ~/nextflow main.nf --input \"path/fastas/\"
		operators: --help calls this help option\n """)
	System.exit(0)
	} else {
		dir= workflow.projectDir.getParent() + "/results/"
		path= params.input + "**/*.fna"
		db = params.db


process loadDB {
}

process identification {
}

process inference {
}

process featureFinder {
}

process visualization {
}
}
