#!/bin/env python

"""Extract information for the QC_visualization out of the PIA output files.
TODO 	nrProteins 	nrProteingroups_unfiltered, number-filtered-protein-groups

DONE number-filtered-peptides, psmZ1-5, psm_missed 0-3 ,number-filtered-psms,
"""

# %% imports
import pandas as pd
import csv
import argparse
import zipfile
import base64
import os
import io

# %% functions
def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pia_peptides", help="peptides.txt from PIA output")
    parser.add_argument("--pia_proteins", help="Proteins.mzID from PIA output")
    parser.add_argument("--pia_PSMs", help="PSM.mzTab from PIA output")
    parser.add_argument("--output", help="path where output should be saved")

    return parser.parse_args()


def run_pia_extraction():
    args = argparse_setup()
    number_proteins, number_ungrouped_proteins = count_nr_Proteins(args.pia_proteins)
    peptide_count = count_nr_filtered_peptides(args.pia_peptides)
    PSM_counts , charge_counts, miss_counts = read_mzTab(args.pia_PSMs)
    dics = [number_proteins, number_ungrouped_proteins, peptide_count, PSM_counts, charge_counts, miss_counts]

    # geschrieben von Dirk (er muss debuggen falls kaputt geht)
    data = {
        key: [value]
        for d in dics
        for key, value in d.items()
    }
    df = pd.DataFrame(data=data)

    zipfile.ZipFile("pia_extractions.zip", mode="w", compresslevel=9).write(args.pia_peptides,
                                                                       compress_type=zipfile.ZIP_DEFLATED,
                                                                       compresslevel=9,
                                                                       arcname=args.pia_peptides.split(os.sep)[-1])
    zipfile.ZipFile("pia_extractions.zip", mode="w", compresslevel=9).write(args.pia_PSMs,
                                                                       compress_type=zipfile.ZIP_DEFLATED,
                                                                       compresslevel=9,
                                                                       arcname=args.pia_PSMs.split(os.sep)[-1])
    zipfile.ZipFile("pia_extractions.zip", mode="w", compresslevel=9).write(args.pia_proteins,
                                                                       compress_type=zipfile.ZIP_DEFLATED,
                                                                       compresslevel=9,
                                                                       arcname=args.pia_proteins.split(os.sep)[-1])
    with open("pia_extractions.zip", "rb") as pia_b:
        pia_str_bs64 = base64.b64encode(pia_b.read())
        df["pia_output.zip"] = pia_str_bs64
        df.to_csv(args.output, index=False)


def count_nr_filtered_peptides(file) -> int:
    # needs the PIA peptides.csv output file and counts PEPTIDE occurences
    # peptides are already filtered for FDR and decoys by the PIA export
    with open(file) as f:
        text = "\n".join([line for line in f if line.startswith("PEPTIDE") or line.startswith("COLS_PEPTIDE") or line.startswith('"COLS_PEPTIDE"')])

    if text:
        peptide_df = pd.read_csv(io.StringIO(text), sep=",")
        
        peptide_count = {}
        peptide_count["number_filtered_peptides"] = peptide_df.shape[0]
    else:
        peptide_count = {}
        peptide_count["number_filtered_peptides"] = 0
    return peptide_count


def count_nr_Proteins(file):
    # read in proteins
    with open(file) as f:
        text = "\n".join([line for line in f if line.startswith("PRH") or line.startswith("PRT")])

    if text:
        protein_df = pd.read_csv(io.StringIO(text), sep="\t")

        # proteins are already FDR and decoy filtered by json parameters
        number_proteins = {}
        number_proteins["number_proteins"] = protein_df.shape[0]

        ungrouped_proteins = list()
        ungrouped_proteins.extend(protein_df['accession'].values.tolist())
        ungrouped_proteins.extend(protein_df[pd.isna(protein_df['ambiguity_members']) == False]['ambiguity_members'].values.tolist())
        number_ungrouped_proteins = {}
        number_ungrouped_proteins["number_ungrouped_proteins"] = len(set(ungrouped_proteins))
    else:
        number_proteins = {}
        number_proteins["number_proteins"] = 0
        number_ungrouped_proteins = {}
        number_ungrouped_proteins["number_ungrouped_proteins"] = 0

    return number_proteins, number_ungrouped_proteins


def read_mzTab(file):
    # get header for PSM FDR Score
    with open(file) as f:
        text = "\n".join([line for line in f if line.startswith("MTD")])

    if text: 
        mzTab_header = pd.read_csv(io.StringIO(text), sep="\t", header=None)
        psm_score_header = mzTab_header[mzTab_header[2].str.contains('MS:1002355')].iloc[-1][1]
        psm_score_header = psm_score_header[4:]     # remove "psm_" at beginning

        ## read in PSMs
        with open(file) as f:
            text = "\n".join([line for line in f if line.startswith("PS")])

        psm_df = pd.read_csv(io.StringIO(text), sep="\t")

        # remove decoys and filter PSM columns
        psm_df = psm_df.loc[psm_df['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
        psm_df = psm_df.loc[:,["PSM_ID", "sequence", "accession", "unique", "retention_time", "charge", "opt_global_missed_cleavages", "modifications", "retention_time", "exp_mass_to_charge", "calc_mass_to_charge", "spectra_ref", psm_score_header]]

        # group the accessions
        gbseries = psm_df.groupby(by=['PSM_ID'])['accession']
        psm_df['accession'] = psm_df["PSM_ID"].map(gbseries.apply(",".join))
        # This uses something incompatible in numpy in some specific versions (e.g. 1.23.5):
        # psm_df['accession'] = psm_df.groupby(by=['PSM_ID'])['accession'].transform(lambda x: ",".join(x))
        psm_df = psm_df.drop_duplicates()

        # filter FDR <= 0.01    # TODO: parameterize?
        psm_df = psm_df.loc[psm_df[psm_score_header] <= 0.01]

        PSM_count = {}
        PSM_count["number_filtered_psms"] = psm_df.shape[0]
        
        nr_psms = psm_df.shape[0] if psm_df.shape[0] != 0 else 1
        charge_counts_above5 = psm_df[psm_df['charge'] > 5]['PSM_ID'].count() / nr_psms
        charge_counts_1      = psm_df[psm_df['charge'] == 1]['PSM_ID'].count() / nr_psms
        charge_counts_2      = psm_df[psm_df['charge'] == 2]['PSM_ID'].count() / nr_psms
        charge_counts_3      = psm_df[psm_df['charge'] == 3]['PSM_ID'].count() / nr_psms
        charge_counts_4      = psm_df[psm_df['charge'] == 4]['PSM_ID'].count() / nr_psms
        charge_counts_5      = psm_df[psm_df['charge'] == 5]['PSM_ID'].count() / nr_psms
        charge_counts = {"psm_charge1": charge_counts_1, "psm_charge2": charge_counts_2, "psm_charge3": charge_counts_3, "psm_charge4": charge_counts_4, "psm_charge5": charge_counts_5, "psm_charge_more": charge_counts_above5}

        miss_count_0    = psm_df[psm_df['opt_global_missed_cleavages'] == 0]['PSM_ID'].count() /nr_psms
        miss_count_1    = psm_df[psm_df['opt_global_missed_cleavages'] == 1]['PSM_ID'].count() /nr_psms
        miss_count_2    = psm_df[psm_df['opt_global_missed_cleavages'] == 2]['PSM_ID'].count() /nr_psms
        miss_count_3    = psm_df[psm_df['opt_global_missed_cleavages'] == 3]['PSM_ID'].count() /nr_psms
        miss_count_more = psm_df[psm_df['opt_global_missed_cleavages'] > 3]['PSM_ID'].count() /nr_psms
        miss_counts = {"psm_missed_0": miss_count_0, "psm_missed_1": miss_count_1, "psm_missed_2": miss_count_2, "psm_missed_3": miss_count_3, "psm_missed_more": miss_count_more}

    else:
        PSM_count = {}
        PSM_count["number_filtered_psms"] = 0
        charge_counts = {"psm_charge1": 0, "psm_charge2": 0, "psm_charge3": 0, "psm_charge4": 0, "psm_charge5": 0, "psm_charge_more": 0}
        miss_counts = {"psm_missed_0": 0, "psm_missed_1": 0, "psm_missed_2": 0, "psm_missed_3": 0, "psm_missed_more": 0}

    return PSM_count, charge_counts, miss_counts


# %% call the script
run_pia_extraction()
