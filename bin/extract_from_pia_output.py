#!/bin/env python

"""Extract information for the QC_visualization out of the PIA output files.

1. Tabelle PSM (aus KNIME):
filtered PSM FDR score < 0.01 
Sequence, Acessions, Modifications, Decoy, Charge, MZ, Deltamass, deltappm, retention time, misscleavages, sourceID, spectrum title, scores, scorenames, 
scoreshorts
Javasnippet in dem Kontaminanten entfernt werden

2. Tabelle 
"""
import pandas as pd
import csv
import io
import argparse
import zipfile
import base64
import os
from lxml import etree as ET
from pathlib import Path

PROTEIN_AMBIGUITY_ELEM: bytes = b'<ProteinAmbiguityGroup'
PROTEIN_HYPOTHESIS_ELEM: bytes = b'<ProteinDetectionHypothesis'
PASS_THRESHOLD_ATTR: bytes = b'passThreshold="'
ID_ATTR: bytes = b'id="'

def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pia_peptides", help="peptides.txt from PIA output")
    parser.add_argument("--pia_proteins", help="Proteins.mzID from PIA output")
    parser.add_argument("--pia_PSMs", help="PSM.mzTab from PIA output")
    parser.add_argument("--output", help="path where output should be saved")

    return parser.parse_args()

def run_pia_extraction():
    args = argparse_setup()
    peptide_count = count_nr_filtered_peptides(args.pia_peptides)
    PSM_counts , charge_counts, miss_counts = read_PSM_mzTab(args.pia_PSMs)
    dics = [peptide_count, PSM_counts, charge_counts, miss_counts]

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
                                                                       arcname=args.proteins.split(os.sep)[-1])
    with open("pia_extractions.zip", "rb") as pia_b:
        pia_str_bs64 = base64.b64encode(pia_b.read())
        df["pia_output.zip"] = pia_str_bs64
        df.to_csv(args.output)


def count_nr_filtered_peptides(file: str) -> int:
    """parses the peptide csv. Important to note: the file alternates between lines for peptides and its PSM. Multiple lines of PSM can follow after a line of peptides"""
    peps_only_csv: str = ""
    psms_only_csv: str = ""
    with Path(file).open('r') as pia_f:
        for line in pia_f:
            if line.startswith("PEPTIDE") or line.startswith('"COLS_PEPTIDE"'):
                peps_only_csv += line
            else:
                psms_only_csv += line
    

    pep_df = pd.read_csv(
        io.StringIO(peps_only_csv),
        sep=","
    )
    psm_df = pd.read_csv(
        io.StringIO(psms_only_csv),
        sep=","
    )
    

  #  print(psm_df[["score names", "scores"]])
   # print(psm_df["score names"].values)
    psm_df_scores = psm_df["scores"].str.split(";", expand=True).add_prefix("score_")
   # print(psm_df_scores)

    psm_df["FDRScore"] = psm_df_scores["score_5"]
    PSM_counts_filtered = psm_df.loc[psm_df["FDRScore"].astype(float) < 0.01]
    #print(psm_df)

   # print(len(pep_df))
    #print(len(PSM_counts_filtered["Sequence"].unique()))
    print(pep_df.columns)
    print(pep_df["score names"].values)
    print(pep_df["best scores"])
    print(len(pep_df["Sequence"].unique()))
    return(len(pep_df["Sequence"].unique()))


    #needs the PIA peptides.csv output file and counts PEPTIDE occurences
    pep_df = pd.read_csv(file, sep=",")
    print(pep_df[["COLS_PEPTIDE", "Sequence"]])
    with open(file) as peps:
        pepsi_count = 0
        pepsi = csv.reader(peps)
        for row in pepsi:
            if(row[0] == "PEPTIDE"):
                pepsi_count += 1
        peptide_count ={}
        peptide_count["number_filtered_peptides"] = pepsi_count
        print("counted peptides: ",peptide_count)
        return peptide_count


def read_PSM_mzTab(file: str):
    # combined FDRScore -> MS:1002356, danach suchen
    with open(file) as tabbi:
        header_line = -1
        for idx, line in enumerate(tabbi):
            if line.startswith("MTD") or line.strip() == "":
                continue
            else:
                header_line = idx - 1
                break
        tabbi.seek(0)
        tab = pd.read_csv(
            tabbi,
            sep="\t",
            header=header_line
        )
        #remove all duplicate PSMs
        tab = tab.drop_duplicates(subset=["PSM_ID"])
        PSM_counts_all = len(tab["PSM_ID"])
        print("all unfiltered PSMs: ", PSM_counts_all)
        # filter for score 
        tab = tab.loc[tab["search_engine_score[1]"] < 0.01]

        """this here is for debugging reasons"""
      #  print(tab.columns)
      #  print(tab["search_engine_score[1]"])
      #  print(tab["search_engine_score[2]"])
      #  print(tab["search_engine_score[3]"])
      #  print(tab["search_engine_score[4]"])
        PSM_counts_filtered = tab.loc[tab["search_engine_score[1]"] < 0.01]

        PSM_counts = len(PSM_counts_filtered["PSM_ID"].unique())
        PSM_count = {}
        PSM_count["number_filtered_pms"] = PSM_counts

        # get charge states
        charge_counts_above5 = 0
        charge_counts_1 = len(tab[tab.charge == 1])
        charge_counts_2 = len(tab[tab.charge == 2])
        charge_counts_3 = len(tab[tab.charge == 3])
        charge_counts_4 = len(tab[tab.charge == 4])
        charge_counts_5 = len(tab[tab.charge >= 5])

        
        #get misscleavages
        miss_count_0 = (len(tab[tab.opt_global_missed_cleavages == 0]) * 100)/ PSM_counts
        miss_count_1 = (len(tab[tab.opt_global_missed_cleavages == 1])* 100)/ PSM_counts
        miss_count_2 = (len(tab[tab.opt_global_missed_cleavages == 2])* 100)/ PSM_counts
        miss_count_3 = (len(tab[tab.opt_global_missed_cleavages == 3])* 100)/ PSM_counts
        miss_count_more = (len(tab[tab.opt_global_missed_cleavages >= 4])* 100)/ PSM_counts
   
        
    

        charge_counts = {"Z1": charge_counts_1, "Z2": charge_counts_2, "Z3": charge_counts_3, "Z4": charge_counts_4, "Z5": charge_counts_5, "Z_more": charge_counts_above5}
        miss_counts = {"missed_0": miss_count_0, "missed_1": miss_count_1, "missed_2": miss_count_2, "missed_3": miss_count_3, "missed_more": miss_count_more}

        print("filtered PSM: ", PSM_counts)
        print(miss_counts)
        return PSM_count, charge_counts, miss_counts

def parse_protein_mzid(xmlfile: Path) -> int:
    """parses mzid based on Index"""

    with xmlfile.open("rb") as f:
        xmlb = f.read()
        num_groups = xmlb.count(PROTEIN_AMBIGUITY_ELEM)
        prots = set()
        next_hypo_offset = 0
        while True:
            # geschrieben von Dirk (er muss debuggen falls kaputt geht)
            # find protein hypo
            start = xmlb.find(PROTEIN_HYPOTHESIS_ELEM, next_hypo_offset)
            if start < 0:
                break
            end = xmlb.find(b">", start)
            if end < 0:
                raise IndexError("Did not find the closing tag")
            next_hypo_offset = end
            hypo = xmlb[start:end]
            # find id attr
            id_start = hypo.find(ID_ATTR)
            id_end = hypo.find(b'"', id_start + len(ID_ATTR))  # start after `id="`
            id = hypo[id_start:id_end]
            if b"DECOY" in id:
                continue
            # find threshold pass
            pass_start = hypo.find(PASS_THRESHOLD_ATTR)
            pass_end = hypo.find(b'"', pass_start + len(PASS_THRESHOLD_ATTR))  # start after `passThreshold="`
            pass_threshold = hypo[pass_start + len(PASS_THRESHOLD_ATTR):pass_end] == b'true'
            if not pass_threshold:
                continue
            prots.add(id.split(b"_")[1])


        print("amount of proteins: ", len(prots))
        print("amount of protein groups :", num_groups)
        return(len(prots))

def parse_protein_mztab(file: str):
        """ not completed, problem with second header not resolved"""
        with open(file) as tabbi:
            header_line = -1
            for idx, line in enumerate(tabbi):
                if line.startswith("MTD") or line.strip() == "":
                    continue
                else:
                    header_line = idx - 1
                    break
            tabbi.seek(0)

            tab = pd.read_csv(
                tabbi,
                sep="\t",
                header=header_line
            )
            print(tab)




        
# OEI37731std_____pia-compilation-piaExport-PSM.mzTab
# FLI19441std_____pia-compilation-piaExport-PSM.mzTab
#QEXI39419std_____pia-compilation-piaExport-PSM.mzTab

read_PSM_mzTab("/home/maike/Programmierprojekte/hackathon5/FLI19441std_____pia-compilation-piaExport-PSM.mzTab")

count_nr_filtered_peptides("/home/maike/Programmierprojekte/hackathon5/QEXI39419std_____pia-compilation-piaExport-peptides.csv")

#QEXI39419std_____pia-compilation-piaExport--proteins.mzid
#parse_protein_mzid(Path("/home/maike/Programmierprojekte/hackathon5/QEXI39419std_____pia-compilation-piaExport--proteins.mzid"))

#parse_protein_mztab("/home/maike/Programmierprojekte/hackathon5/FLI19441std_____pia-compilation-piaExport--proteins.mzTab")