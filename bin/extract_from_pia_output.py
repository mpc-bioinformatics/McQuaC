#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import io
import h5py
from typing import List, Any


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pia_peptides", help="peptides.txt from PIA output")
    parser.add_argument("--pia_proteins", help="Proteins.mzTab from PIA output")
    parser.add_argument("--pia_PSMs", help="PSM.mzTab from PIA output")
    parser.add_argument("--out_hdf5", help="Output HDF5 with statistics")
    parser.add_argument("--store_all_infos", help="Store all PSMs, peptides and proteins in the HDF5 file", default=False, type=bool)

    return parser.parse_args()


def add_entry_to_hdf5(
    f, qc_acc: str, qc_short_name: str, qc_name: str, qc_description: str, 
    value, value_shape: tuple, value_type: str, 
    unit_accession: str=None, unit_name: str=None,
    ): 
    """ Adds an entry into the hdf5 file """
    key = "|".join([qc_acc, qc_short_name])  # ACCESSION|SHORT_DESC
    if value_type in ("str", h5py.string_dtype()):
        ds = f.create_dataset(key, shape=value_shape, dtype=h5py.string_dtype(), compression="gzip")
        ds[:] = value
    else:
        f.create_dataset(key, value_shape, dtype=value_type, compression="gzip")
        f[key].write_direct(np.array(value, dtype=value_type))
    
    f[key].attrs["qc_short_name"] = qc_short_name
    f[key].attrs["qc_name"] = qc_name
    f[key].attrs["qc_description"] = qc_description
    f[key].attrs["unit_accession"] = unit_accession
    f[key].attrs["unit_name"] = unit_name


def add_table_to_hdf5(
    f, qc_acc: str, qc_short_name: str, qc_name: str, qc_description: str, 
    column_names: List[str], column_data: List[List[Any]], column_types: List[str]
    ): 
    """Adds a table in groups"""

    key = "|".join([qc_acc, qc_short_name])  # ACCESSION|SHORT_DESC
    table_group = f.create_group(key)
    table_group.attrs["qc_short_name"] = qc_short_name
    table_group.attrs["qc_name"] = qc_name
    table_group.attrs["qc_description"] = qc_description
    table_group.attrs["column_order"] = "|".join(column_names)
    
    for n, d, t in zip (column_names, column_data, column_types): 
        if t in ("str", h5py.string_dtype()):
            ds = table_group.create_dataset(n, shape=len(d), dtype=h5py.string_dtype(), compression="gzip")
            ds[:] = d
        else:
            table_group.create_dataset(n, (len(d),), dtype=t, compression="gzip")
            table_group[n].write_direct(np.array(d, dtype=t))


def run_pia_extraction():
    args = argparse_setup()

    # Open HDF5 file in write mode
    with h5py.File(args.out_hdf5, 'w') as out_h5:
        parse_protein_infos(args.pia_proteins, out_h5)
        parse_peptide_infos(args.pia_peptides, out_h5)
        parse_psm_infos(args.pia_PSMs, out_h5, args.store_all_infos)


def parse_peptide_infos(pia_peptide_csv: str, out_hdf5: h5py.File) -> int:
    # needs the PIA peptides.csv output file and counts PEPTIDE occurences
    # peptides are already filtered for FDR and decoys by the PIA export
    with open(pia_peptide_csv) as f:
        text = "\n".join([line for line in f if line.startswith("PEPTIDE") or line.startswith("COLS_PEPTIDE") or line.startswith('"COLS_PEPTIDE"')])

    if text:
        peptide_df = pd.read_csv(io.StringIO(text), sep=",")
        number_filtered_peptides = peptide_df.shape[0]
    else:
        number_filtered_peptides = 0
    
    add_entry_to_hdf5(out_hdf5,
                      "MS:1003250", "number_peptides", "count of identified peptidoforms", "The number of peptidoforms that pass the threshold to be considered identified with sufficient confidence.", 
                      number_filtered_peptides, (1,), "int32", 
                      "UO:0000189", "count unit")


def parse_protein_infos(pia_proteins_mztab: str, out_hdf5: h5py.File):
    # read in proteins
    with open(pia_proteins_mztab) as f:
        text = "\n".join([line for line in f if line.startswith("PRH") or line.startswith("PRT")])

    if text:
        protein_df = pd.read_csv(io.StringIO(text), sep="\t")

        # proteins are already FDR and decoy filtered by PIA settings
        number_proteins = protein_df.shape[0]

        ungrouped_proteins = list()
        ungrouped_proteins.extend(protein_df['accession'].values.tolist())
        ungrouped_proteins.extend(protein_df[pd.isna(protein_df['ambiguity_members']) == False]['ambiguity_members'].values.tolist())
        number_ungrouped_proteins = len(set(ungrouped_proteins))
    else:
        number_proteins = 0
        number_ungrouped_proteins = 0
    
    add_entry_to_hdf5(out_hdf5,
                      "MS:1003327", "number_proteins", "number of identified protein groups", "The number of protein groups that pass the threshold to be considered identified with sufficient confidence.", 
                      number_proteins, (1,), "int32", 
                      "UO:0000189", "count unit")
    
    add_entry_to_hdf5(out_hdf5,
                      "Local:08", "number_ungrouped_proteins", "number of ungrouped identified protein groups", "TODO", 
                      number_ungrouped_proteins, (1,), "int32", 
                      "UO:0000189", "count unit")


def parse_psm_infos(pia_psm_mztab: str, out_hdf5: h5py.File, store_all_infos: bool = False):
    # get header for PSM FDR Score
    with open(pia_psm_mztab) as f:
        text = "\n".join([line for line in f if line.startswith("MTD")])

    if text: 
        mzTab_header = pd.read_csv(io.StringIO(text), sep="\t", header=None)
        psm_score_header = mzTab_header[mzTab_header[2].str.contains('MS:1002355')].iloc[-1][1] # this is the PSM-level FDR score
        psm_score_header = psm_score_header[4:]     # remove "psm_" at beginning

        ## read in PSMs
        with open(pia_psm_mztab) as f:
            text = "\n".join([line for line in f if line.startswith("PS")])

        psm_df = pd.read_csv(io.StringIO(text), sep="\t")

        # remove decoys and filter for important PSM columns
        psm_df = psm_df.loc[psm_df['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
        psm_df = psm_df.loc[:,["PSM_ID", "sequence", "accession", "unique", "retention_time", "charge", "opt_global_missed_cleavages", "modifications", "exp_mass_to_charge", "calc_mass_to_charge", "spectra_ref", psm_score_header]]
        
        # group the accessions
        gbseries = psm_df.groupby(by=['PSM_ID'])['accession']
        psm_df['accession'] = psm_df["PSM_ID"].map(gbseries.apply(",".join))
        psm_df = psm_df.drop_duplicates()

        unfiltered_psms = psm_df.copy()
        unfiltered_psms.rename(columns={psm_score_header: 'psm_level_fdrscore'}, inplace=True)

        # filter FDR <= 0.01    # TODO: parameterize?
        psm_df = psm_df.loc[psm_df[psm_score_header] <= 0.01]

        # Calculate ppm error
        exp_calc_diff = psm_df["exp_mass_to_charge"] - psm_df["calc_mass_to_charge"]
        exp_calc_diff_removed_isotopes = exp_calc_diff - (exp_calc_diff.round()).astype(int) # Remove Isotopes, since calc_mass expects none
        ppm_error_df = (exp_calc_diff_removed_isotopes * 1000000) / psm_df["calc_mass_to_charge"]
        ppm_error = ppm_error_df.to_list()

        # now get all counts as QC metrics
        PSM_count = psm_df.shape[0]
        
        nr_psms = psm_df.shape[0] if psm_df.shape[0] != 0 else 1
        charge_fraction_1    = psm_df[psm_df['charge'] == 1]['PSM_ID'].count() / nr_psms
        charge_fraction_2    = psm_df[psm_df['charge'] == 2]['PSM_ID'].count() / nr_psms
        charge_fraction_3    = psm_df[psm_df['charge'] == 3]['PSM_ID'].count() / nr_psms
        charge_fraction_4    = psm_df[psm_df['charge'] == 4]['PSM_ID'].count() / nr_psms
        charge_fraction_5    = psm_df[psm_df['charge'] == 5]['PSM_ID'].count() / nr_psms
        charge_fraction_more = psm_df[psm_df['charge'] > 5]['PSM_ID'].count() / nr_psms
        charge_fractions = [[charge_fraction_1], [charge_fraction_2], [charge_fraction_3], [charge_fraction_4], [charge_fraction_5], [charge_fraction_more]]

        miss_count_0    = psm_df[psm_df['opt_global_missed_cleavages'] == 0]['PSM_ID'].count()
        miss_count_1    = psm_df[psm_df['opt_global_missed_cleavages'] == 1]['PSM_ID'].count()
        miss_count_2    = psm_df[psm_df['opt_global_missed_cleavages'] == 2]['PSM_ID'].count()
        miss_count_3    = psm_df[psm_df['opt_global_missed_cleavages'] == 3]['PSM_ID'].count()
        miss_count_more = psm_df[psm_df['opt_global_missed_cleavages'] > 3]['PSM_ID'].count()
        miss_counts = [[miss_count_0], [miss_count_1], [miss_count_2], [miss_count_3], [miss_count_more]]

    else:
        PSM_count = 0
        charge_fractions = [[0], [0], [0], [0], [0], [0]]
        miss_counts = {"psm_missed_0": 0, "psm_missed_1": 0, "psm_missed_2": 0, "psm_missed_3": 0, "psm_missed_more": 0}
        ppm_error = [np.nan]

    add_entry_to_hdf5(out_hdf5,
                      "MS:1003251", "number_of_filtered_psms", "count of identified spectra", "The number of spectra that pass the threshold to be considered identified with sufficient confidence.", 
                      PSM_count, (1,), "int32", 
                      "UO:0000189", "count unit")
    
    add_table_to_hdf5(out_hdf5,
                      "Local:09", "psm_charge_fractions", "TODO", "TODO", 
                      ["1", "2", "3", "4", "5", "6 or more"], charge_fractions, ["float64", "float64", "float64", "float64", "float64", "float64"]
    )

    add_table_to_hdf5(out_hdf5,
                      "MS:4000180", "psm_missed_counts", "table of missed cleavage counts", "The number of identified peptides with corresponding number of missed cleavages after user-defined acceptance criteria are applied. The number of missed cleavages per peptide is given in the 'number of missed cleavages' column, the respective count of such peptides identified in the 'Number of Occurrences' column. The highest 'missed cleavages' row is to be interpreted as that number of missed cleavages or higher.", 
                      ["0", "1", "2", "3", "4 or more"], miss_counts, ["uint32", "uint32", "uint32", "uint32", "uint32"]
    )

    add_entry_to_hdf5(out_hdf5,
                      "LOCAL:05", "filtered_psms_ppm_error", "deviations in PPM for each PSM", "TODO", 
                      ppm_error, (len(ppm_error),), "float64", 
                      "UO:0000169", "parts per million")
    
    if store_all_infos and (unfiltered_psms is not None):
        # TODO: fix in PIA
        unfiltered_psms["modifications"] = unfiltered_psms["modifications"].astype("str")

        column_names = unfiltered_psms.columns.tolist()
        psms_data = [unfiltered_psms[x].astype("str") for x in column_names]
        column_types = unfiltered_psms.dtypes.astype("str").tolist()

        add_table_to_hdf5(out_hdf5,
                        "identified_psms_table", "identified_psms_table", "table of all unfiltered identified PSMs", "All identified PSM as reported by PIA", 
                        column_names, psms_data, column_types
        )

if __name__ == "__main__":
    run_pia_extraction()
