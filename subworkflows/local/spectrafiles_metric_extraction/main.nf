include { MZMLMETRICSEXTRACTION } from '../../../modules/local/mzmlmetricsextraction/main'

workflow SPECTRAFILES_METRIC_EXTRACTION {

    take:
    mzmls           // channel: [ val(meta), mzml ]
    spectra_files   // channel: [ val(meta), spectra_file ]

    main:

    // first extract data from mzML files
    MZMLMETRICSEXTRACTION(
        mzmls
    )

    emit:
    mzml_metrics      = MZMLMETRICSEXTRACTION.out.hdf5
}
