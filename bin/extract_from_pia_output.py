#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import zipfile
import base64
import os
import io
import h5py


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pia_peptides", help="peptides.txt from PIA output")
    parser.add_argument("--pia_proteins", help="Proteins.mzID from PIA output")
    parser.add_argument("--pia_PSMs", help="PSM.mzTab from PIA output")
    parser.add_argument("--out_hdf5", help="Output HDF5 with statistics")

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



def add_table_in_hdf5(
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
        count_nr_Proteins(args.pia_proteins, out_h5)
        count_nr_filtered_peptides(args.pia_peptides, out_h5)
        read_mzTab(args.pia_PSMs, out_h5)

        # Save the zip files also in HDF5
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
            pia_str_bs64 = base64.b64encode(pia_b.read()).decode("utf-8")
        out_h5.create_dataset("pia_output.zip", data=pia_str_bs64)
        out_h5["pia_output.zip"].attrs["unit"] = "zip File"
        out_h5["pia_output.zip"].attrs["Description"] = "A zip container containing the used PSM-mzTAB, Peptide-CSV and Proteins-mzTAB file from the PIA-Output. This zip file was base64 encoded to be represented as a string"


def count_nr_filtered_peptides(pia_peptide, out_hdf5) -> int:
    # needs the PIA peptides.csv output file and counts PEPTIDE occurences
    # peptides are already filtered for FDR and decoys by the PIA export
    with open(pia_peptide) as f:
        text = "\n".join([line for line in f if line.startswith("PEPTIDE") or line.startswith("COLS_PEPTIDE") or line.startswith('"COLS_PEPTIDE"')])

    if text:
        peptide_df = pd.read_csv(io.StringIO(text), sep=",")
        number_filtered_peptides = peptide_df.shape[0]
    else:
        number_filtered_peptides = 0

    add_entry_to_hdf5(
        out_hdf5, "number_filtered_peptides", number_filtered_peptides, (1,), "int32", "none", 
        description="Peptides filtered by an FDR-Score of 0.01"
    )


def count_nr_Proteins(pia_prots, out_hdf5):
    # read in proteins
    with open(pia_prots) as f:
        text = "\n".join([line for line in f if line.startswith("PRH") or line.startswith("PRT")])

    if text:
        protein_df = pd.read_csv(io.StringIO(text), sep="\t")

        # proteins are already FDR and decoy filtered by json parameters
        number_proteins = protein_df.shape[0]

        ungrouped_proteins = list()
        ungrouped_proteins.extend(protein_df['accession'].values.tolist())
        ungrouped_proteins.extend(protein_df[pd.isna(protein_df['ambiguity_members']) == False]['ambiguity_members'].values.tolist())
        number_ungrouped_proteins = len(set(ungrouped_proteins))
    else:
        number_proteins = 0
        number_ungrouped_proteins = 0

    add_entry_to_hdf5(
        out_hdf5, "number_proteins", number_proteins, (1,), "int32", "none", 
        description="Number of protein groups."
    )
    add_entry_to_hdf5(
        out_hdf5, "number_ungrouped_proteins", number_ungrouped_proteins, (1,), "int32", "none", 
        description="Total number of protein accessions within the protein groups."
    )


def read_mzTab(pia_mzTAB, out_hdf5):
    # get header for PSM FDR Score
    with open(pia_mzTAB) as f:
        text = "\n".join([line for line in f if line.startswith("MTD")])

    if text: 
        mzTab_header = pd.read_csv(io.StringIO(text), sep="\t", header=None)
        psm_score_header = mzTab_header[mzTab_header[2].str.contains('MS:1002355')].iloc[-1][1]
        psm_score_header = psm_score_header[4:]     # remove "psm_" at beginning

        ## read in PSMs
        with open(pia_mzTAB) as f:
            text = "\n".join([line for line in f if line.startswith("PS")])

        psm_df = pd.read_csv(io.StringIO(text), sep="\t")

        # remove decoys and filter PSM columns
        psm_df = psm_df.loc[psm_df['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
        psm_df = psm_df.loc[:,["PSM_ID", "sequence", "accession", "unique", "retention_time", "charge", "opt_global_missed_cleavages", "modifications", "retention_time", "exp_mass_to_charge", "calc_mass_to_charge", "spectra_ref", psm_score_header]]
        

        ### Calculate ppm error
        exp_calc_diff = (psm_df["exp_mass_to_charge"] - psm_df["calc_mass_to_charge"])
        exp_calc_diff_removed_isotopes = (exp_calc_diff - (exp_calc_diff.round()).astype(int)) # Remove Isotopes, since calc_mass expects none
        ppm_error_df = ((exp_calc_diff_removed_isotopes * 1000000) / psm_df["calc_mass_to_charge"])
        ppm_error = {"filtered_psms_ppm_error": ppm_error_df.to_list()}
        # Get the ppm error to the theoretical masses

        (psm_df["calc_mass_to_charge"] - exp_calc_diff_removed_isotopes) # measured
        psm_df["calc_mass_to_charge"] # expected

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
        PSM_count = {"number_filtered_psms": 0}
        charge_counts = {"psm_charge1": 0, "psm_charge2": 0, "psm_charge3": 0, "psm_charge4": 0, "psm_charge5": 0, "psm_charge_more": 0}
        miss_counts = {"psm_missed_0": 0, "psm_missed_1": 0, "psm_missed_2": 0, "psm_missed_3": 0, "psm_missed_more": 0}
        ppm_error = {"filtered_psms_ppm_error": [np.nan]}

    # concat all entries into a single dict
    dicts = [PSM_count, charge_counts, miss_counts]
    data = {
        key: [value]
        for d in dicts
        for key, value in d.items()
    }
    
    # Keys are added in order in dict, therefore we could leave out the order, but set this order statically to be save just in case!!!
    keys = ["number_filtered_psms","psm_charge1","psm_charge2","psm_charge3","psm_charge4","psm_charge5","psm_charge_more",
        "psm_missed_0","psm_missed_1","psm_missed_2","psm_missed_3","psm_missed_more"]
    descriptions = [
        "PSMs filtered by an FDR-Score of 0.01",
        "proportion of identified MS2 spectra with charge 1",
        "proportion of identified MS2 spectra with charge 2",
        "proportion of identified MS2 spectra with charge 3",
        "proportion of identified MS2 spectra with charge 4",
        "proportion of identified MS2 spectra with charge 5",
        "proportion of identified MS2 spectra with charge 6 or more",
        "proportion of PSMs with a sequence containing 0 missed cleavages",
        "proportion of PSMs with a sequence containing 1 missed cleavages",
        "proportion of PSMs with a sequence containing 2 missed cleavages",
        "proportion of PSMs with a sequence containing 3 missed cleavages",
        "proportion of PSMs with a sequence containing 4 ore more missed cleavages"
    ]

    for k, d in zip(keys, descriptions):
        add_entry_to_hdf5(
            out_hdf5, k, data[k], (1,), "float64", "none", 
            description=d
        )

    if len(ppm_error["filtered_psms_ppm_error"]) > 0:
        add_entry_to_hdf5(
            out_hdf5, "filtered_psms_ppm_error", ppm_error["filtered_psms_ppm_error"], (len(ppm_error["filtered_psms_ppm_error"]),), "float64", "ppm", 
            description="PPM-Error from the calculated (theoretical) to the experimental (measured) precursor. We calculated here 'theoretical - calc'. "
        )
    else:
        add_entry_to_hdf5(
            out_hdf5, "filtered_psms_ppm_error", np.nan, (1,), "float64", "ppm", 
            description="PPM-Error from the calculated (theoretical) to the experimental (measured) precursor. We calculated here 'theoretical - calc'. "
        )

if __name__ == "__main__":
    run_pia_extraction()
