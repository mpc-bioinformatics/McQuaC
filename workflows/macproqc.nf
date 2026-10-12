/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_macproqc_pipeline'

include { PREPARE_SPECTRA } from '../subworkflows/local/prepare_spectra'
include { IDENT_DDA } from '../subworkflows/local/ident_dda'
include { SPECTRAFILES_METRIC_EXTRACTION } from '../subworkflows/local/spectrafiles_metric_extraction'

include { COMBINEHDF5 } from '../modules/local/combinehdf5/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MACPROQC {

    take:
    ch_samplesheet              // channel: samplesheet read in from --input
    outdir
    fasta
    skip_decoy_generation
    comet_config_template       // string: path to comet config template, or null to use created default config
    label_modifications         //
    pia_fdr_threshold           // number: pia fdr threshold (default 0.01)
    pia_prefilter_threshold     // number: pia prefilter threshold (the FDR pre-filter or 0 / empty / null to skip prefiltering)

    main:

    def ch_versions = channel.empty()   // this channel needs to be filled, if the module/subworkflow return any versions NOT in the topic channel, but in the versions.yml file

    // prepare the spectra files (strip the fasta_file before passing to PREPARE_SPECTRA)
    PREPARE_SPECTRA(
        ch_samplesheet.map { meta, spectrum_file -> [meta, spectrum_file] }
    )

    // enable "label search" if the user has specified any label modifications (see also modules.config)
    def search_label_modifications = label_modifications?.trim() ? true : false

    // perform DDA identification and extract the respective metrics
    IDENT_DDA(
        fasta,
        skip_decoy_generation,
        comet_config_template,
        search_label_modifications,
        PREPARE_SPECTRA.out.mzmls,
        pia_fdr_threshold,
        pia_prefilter_threshold
    )

    // extract metrics from the spectra files (mzML, Thermo .raw, Bruker .d)
    SPECTRAFILES_METRIC_EXTRACTION(
        PREPARE_SPECTRA.out.mzmls,
        PREPARE_SPECTRA.out.uncompressed
    )

    // collect all metric HDF5 files per run: the identification metrics carry "<run>-unlabelled"
    // as meta.id (and the run ID in meta.mzml_id), so join them on the run ID and keep the run's meta
    ch_run_metrics = SPECTRAFILES_METRIC_EXTRACTION.out.mzml_metrics
        .map { meta, hdf5 -> [ meta.id, meta, hdf5 ] }
        .join(
            IDENT_DDA.out.idmetrics.map { meta, hdf5 -> [ meta.mzml_id, hdf5 ] },
            remainder: true     // keep runs without identification metrics
        )
        .map { _run_id, meta, mzml_hdf5, id_hdf5 -> [ meta, [ mzml_hdf5, id_hdf5 ].findAll { hdf5 -> hdf5 != null } ] }

    // combine the metrics of each run into one HDF5 file
    COMBINEHDF5(
        ch_run_metrics
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'macproqc_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
