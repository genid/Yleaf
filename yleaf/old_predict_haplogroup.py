#!/usr/bin/env python

"""
Script for Haplogroup prediction

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Diego Montiel Gonzalez
"""

import pandas as pd
import numpy as np
import collections
import operator
import os
import argparse
from argparse import ArgumentParser

from yleaf import yleaf_constants


def get_arguments():
    parser = ArgumentParser(description="Erasmus MC: Genetic Identification\n Y-Haplogroup Prediction")

    parser.add_argument("-input", "--Input",
                        dest="Input", required=True, type=file_exists,
                        help="Output file or path produced from Yleaf", metavar="FILE")

    parser.add_argument("-out", "--Outfile",
                        dest="Outputfile", required=True,
                        help="Output file name", metavar="FILE")

    args = parser.parse_args()
    return args


def file_exists(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError(f"{x} does not exist")
    return x


def check_if_folder(path, ext):
    list_files = []
    if os.path.isdir(path):
        dirpath = os.walk(path)
        for dirpath, dirnames, filenames in dirpath:
            for filename in [f for f in filenames if f.endswith(ext)]:
                files = os.path.join(dirpath, filename)
                list_files.append(files)
        return list_files
    else:
        return [path]


def get_hg_root(hg):
    """
    Choose the haplogroup based on the highest count and branch depth (E.x. J2a21 = J)
    """
    list_hg = []
    init_hg = "NA"
    for i in hg:
        list_hg.append(i)
    collections_hg = collections.Counter(list_hg)
    try:
        hg_dict = {}
        for c in collections_hg:
            if c[0] not in hg_dict:
                hg_dict[c[0]] = collections_hg[c]
            else:
                hg_dict[c[0]] += collections_hg[c]
        init_hg = max(iter(hg_dict.items()), key=operator.itemgetter(1))[0]
        return init_hg
    except:
        return init_hg


def get_intermediate_branch(init_hg, path_hg_prediction_tables):
    tmp_init_hg = init_hg + "_int.txt"
    hg_intermediate_file = path_hg_prediction_tables + "/" + tmp_init_hg
    try:
        df_intermediate = pd.read_csv(hg_intermediate_file, header=None, sep="\t", engine='python')
        return df_intermediate
    except:
        return pd.DataFrame()


def calc_score_one(df_intermediate, df_haplogroup):
    """
    QC.1: Calculates the estimate of correct states from intermediate
    table and the total of intermediates found
    """
    total = 0
    correct_state = 0
    for i in df_intermediate.values:
        tmp = df_haplogroup.loc[df_haplogroup["haplogroup"] == i[0]]
        if not tmp.empty:
            if "/" in i[1]:
                correct_state += len(tmp)
                total += len(tmp)
            else:
                correct_state += np.sum(i[1] == tmp["state"])
                total += len(tmp)
    try:
        qc_one = round((correct_state / total), 3)
    except ZeroDivisionError:
        qc_one = 0.0
    return qc_one


def get_putative_hg_list(init_hg, hg, df_haplogroup_trimmed, df_haplogroup_all):
    """
    Removes all haplogroup but the main one
    Check if the preffix of all main haplogroup with D state by allowing one haplogroup that does not match
    to the main haplogroup until you get back to the origin.
    If mismatch > 1 skip choose as main haplogroup the following in the largest length and continue.
    In case haplogroup contains ~ symbol at the end ignore it temporally for preffix comparison but
    if this is the main haplogroup store the special character at the end.
    Ancestral and Derive state
    QC.2
    Check if the same name of the main haplogroup appears as an Ancestral State and
    save the number of count and calculate QC2
    """
    dict_hg = {}
    list_hg_remove = []
    hg_threshold = 0.8
    list_main_hg = sorted(list(set(hg)), reverse=False)

    for putative_hg in list_main_hg:
        # print(putative_hg)
        total_qctwo = len(df_haplogroup_all.loc[df_haplogroup_all["haplogroup"] == putative_hg])
        ahg = np.sum("A" == df_haplogroup_all.loc[df_haplogroup_all["haplogroup"] == putative_hg]["state"])
        try:
            qc_two = round((total_qctwo - ahg) / total_qctwo, 3)
        except ZeroDivisionError:
            qc_two = 0.0
        if qc_two >= hg_threshold:
            qc_three = calc_score_three(df_haplogroup_trimmed, putative_hg, init_hg)
            dict_hg[putative_hg] = [qc_two, qc_three]

    if not bool(dict_hg):
        return dict_hg
    # in case of a A haplogroup
    if init_hg == 'A':
        key = max(dict_hg)
        dict_hg.update({key: [dict_hg[key][0], 1]})
    else:
        # removes hg with lower qc
        for key in dict_hg:
            if dict_hg[key][0] < hg_threshold or dict_hg[key][1] < hg_threshold:
                list_hg_remove.append(key)
        if len(list_hg_remove) > 0:
            for key in list_hg_remove:
                if len(key) > 1:
                    dict_hg.pop(key, None)
    return dict_hg


def get_putative_hg(dict_hg):
    list_putative_hg = list(dict_hg.keys())
    list_putative_hg.sort(reverse=True)
    putative_hg = "NA"
    qc_two = 0.0
    # Check for same pattern in the derive state
    for i in range(len(list_putative_hg) - 1):
        if list_putative_hg[i + 1] not in list_putative_hg[i]:
            dict_hg.pop(list_putative_hg[i], None)
    if len(list_putative_hg) > 0:
        putative_hg = max(dict_hg.keys(), key=len)
        qc_two = dict_hg[putative_hg]

    return putative_hg, qc_two


def calc_score_three(df_haplogroup, putative_hg, init_hg):
    """
    QC.3
    Show both Ancestral and Derived states from the main haplogroup and check the preffix
    of each haplogroup that is before the main one until you reach the root. If there are
    some Ancestral states which follows the pattern/preffix from the main haplogroup keep
    the count of how many of these appear and this will give the total count. Substract
    the ones found from the corrected and this will give the QC.3 score.
    """
    list_main_hg_all = []
    list_hg_all = df_haplogroup["haplogroup"].values
    for i in list_hg_all:
        if i.startswith(init_hg):
            list_main_hg_all.append(i)
    list_main_hg_all = sorted(list(set(list_main_hg_all)), reverse=True)
    df_main_hg_all = df_haplogroup.loc[df_haplogroup["haplogroup"].isin(list_main_hg_all)]
    df_main_hg_all = df_main_hg_all[["haplogroup", "state", "marker_name"]]
    df_main_hg_all = df_main_hg_all.sort_values(by="haplogroup", ascending=False).values
    a_match = 0
    total_match = 0
    for putative_hg in [putative_hg]:
        for i in df_main_hg_all:
            tmp_hg = i[0]
            if tmp_hg in putative_hg and tmp_hg != putative_hg:
                total_match += 1
                if i[1] == "A":
                    a_match += 1
    try:
        qc_three = round((total_match - a_match) / total_match, 3)
    except ZeroDivisionError:
        qc_three = 0.0
    return qc_three


def get_putative_ancenstral_hg(df_haplogroup, putative_hg):
    """
    Haplogroup and marker name
    if main haplogroup is a preffix from a higher resolution Ancestral state haplogroup
    and report this haplogroup and marker name (Ex. HG: J21a Marker: L123)
    Could be that that there are more than one contain as a preffix from an ancestral state haplogroup.
    Should report all of them only if there are different haplogroup name with the resolution
    """
    putative_ancestral_hg = []
    putative_hg = putative_hg.replace("~", "")
    df_putative_ancestral_hg = df_haplogroup.copy()
    df_putative_ancestral_hg = df_putative_ancestral_hg[df_putative_ancestral_hg.haplogroup.str.startswith(putative_hg)]
    df_putative_ancestral_hg['haplogroup'] = df_putative_ancestral_hg['haplogroup'].str.replace('~', '')
    df_putative_ancestral_hg = df_putative_ancestral_hg[~df_putative_ancestral_hg.haplogroup.isin([putative_hg])]
    df_putative_ancestral_hg = df_putative_ancestral_hg[df_putative_ancestral_hg.state == "A"]
    df_putative_ancestral_hg = df_putative_ancestral_hg.sort_values(by=['haplogroup'])

    for i in df_putative_ancestral_hg.index:
        if len(putative_ancestral_hg) == 0:
            putative_ancestral_hg.append(df_putative_ancestral_hg.loc[i])
        else:
            if putative_ancestral_hg[-1][3] not in df_putative_ancestral_hg.loc[i][3]:
                putative_ancestral_hg.append(df_putative_ancestral_hg.loc[i])
    if len(putative_ancestral_hg) > 0:
        putative_ancestral_hg = pd.DataFrame(putative_ancestral_hg)
    return putative_ancestral_hg


def process_keys(keys):
    tmp_hg = []
    for i in range(len(keys)):
        if "~" in keys[i]:
            tmp_hg.append(keys[i])
    for i in tmp_hg:
        if i in keys:
            index = keys.index(i)
            keys.pop(index)
    return keys, tmp_hg


def process_log(log_file):
    log_file += "info"
    total_reads = "NA"
    valid_markers = "NA"

    try:
        df_log = pd.read_csv(log_file, sep=":", header=None, index_col=0)
        total_reads = df_log.loc["Total of mapped reads"][1]
        valid_markers = df_log.loc["Markers with haplogroup information"][1]
    except FileNotFoundError:
        print("Warning: log file not found!")
    return total_reads, valid_markers


def main():
    print("\tY-Haplogroup Prediction")

    args = get_arguments()

    path_samples = args.Input  # .out files are collected
    samples = check_if_folder(path_samples, '.out')
    out_file = args.Outputfile
    hg_intermediate = str(yleaf_constants.DATA_FOLDER / yleaf_constants.HG_PREDICTION_FOLDER)
    intermediate_tree_table = hg_intermediate + "/Intermediates.txt"
    h_flag = True
    log_output = []
    for sample_name in samples:
        putative_hg = "NA"
        out_name = str(sample_name.split("/")[-1])
        out_name = out_name.split(".")[0]

        total_reads, valid_markers = process_log(sample_name[:-3])

        df_intermediate = pd.read_csv(intermediate_tree_table, header=None, engine='python')
        intermediates = df_intermediate[0].values

        df_haplogroup_all = pd.read_csv(sample_name, sep="\t", engine='python')
        df_haplogroup_all = df_haplogroup_all.sort_values(by=['haplogroup'])

        df_haplogroup_trimmed = df_haplogroup_all.copy()

        df_derived = df_haplogroup_all.copy()
        df_derived = df_derived[df_derived["state"] == "D"]

        df_haplogroup_trimmed['haplogroup'] = df_haplogroup_trimmed['haplogroup'].str.replace('~', '')
        # instance with only D state
        df_tmp = df_derived
        for hg in intermediates:
            # Removes intermediate branches
            df_tmp = df_tmp.drop(df_tmp[df_tmp.haplogroup == hg].index)
        hg = df_tmp["haplogroup"].values

        # pre filter based on the most prevalent and deepest haplogroup in derived data
        init_hg = get_hg_root(hg)

        df_intermediate = get_intermediate_branch(init_hg, hg_intermediate)
        qc_one = calc_score_one(df_intermediate, df_haplogroup_trimmed)

        df_haplogroup_trimmed = df_haplogroup_trimmed[~df_haplogroup_trimmed.haplogroup.isin(intermediates)]

        df_derived = df_derived[~df_derived.haplogroup.isin(intermediates)]
        hg = df_derived[(df_derived.haplogroup.str.startswith(init_hg))].haplogroup.values

        dict_hg = get_putative_hg_list(init_hg, hg, df_haplogroup_trimmed, df_haplogroup_all)
        keys = sorted(dict_hg.keys(), reverse=True, key=lambda x: x.replace("~", ""))  # apparently tilde has a
        # higher ASCII value, resulting in a wrong order

        mismatches = []
        t = 2  # max mismatch for preffix
        # look for the preffix from bottom to the root of the tree
        for k in range(len(keys)):
            mismatch = 0
            for j in range(k + 1, len(keys)):
                if keys[j].replace("~", "") not in keys[k].replace("~", ""):
                    mismatch += 1

            if mismatch < t:
                putative_hg = keys[k]
                qc_two = dict_hg[keys[k]][0]
                qc_three = dict_hg[keys[k]][1]
                break
            mismatches.append(mismatch)

        hg_list = []
        flag = False
        for i in keys:
            if putative_hg.replace("~", "") in i:
                hg_list.append(i)
                flag = True
        if flag:
            putative_hg = max(hg_list, key=len)

        putative_ancestral_hg = get_putative_ancenstral_hg(df_haplogroup_all, putative_hg)

        # Output
        header = "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3"
        marker_name = df_haplogroup_all.loc[df_haplogroup_all["haplogroup"] == putative_hg]["marker_name"].values
        if putative_hg == "NA":
            out_hg = "NA"
            output = "{}\tNA\tNA\t{}\t{}\t0\t0\t0\t0".format(out_name, total_reads, valid_markers)
            log_output.append(out_name)
        else:
            if len(marker_name) > 1:
                out_hg = putative_hg[0] + "-" + marker_name[0] + "/etc"
            elif len(marker_name) == 1:
                out_hg = putative_hg[0] + "-" + marker_name[0]
            if len(putative_ancestral_hg) > 0:
                out_hg += "*(x"
                for i in putative_ancestral_hg.index:
                    out_hg += putative_ancestral_hg.loc[i]["marker_name"] + ","
                out_hg += ")"
                out_hg = list(out_hg)
                del out_hg[-2]
                out_hg = "".join(out_hg)
            qc_score = round((qc_one * qc_two * qc_three), 3)
            if qc_score >= 0.7:
                output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(out_name, putative_hg, out_hg, total_reads,
                                                                     valid_markers, qc_score, qc_one, qc_two, qc_three)
            else:
                log_output.append(out_name)
                output = "{}\tNA\tNA\t{}\t{}\t{}\t{}\t{}\t{}".format(out_name, total_reads, valid_markers, qc_score,
                                                                     qc_one, qc_two, qc_three)

        with open(out_file, "a") as w_file:
            if h_flag:
                h_flag = False
                w_file.write(header)
            w_file.write("\n")
            w_file.write(output)

    if len(log_output) > 0:
        print("Warning: Following sample(s) showed discrepancies, please check output(s) manually: ")
        print("\n".join(log_output))

    print("--- Yleaf 'Y-Haplogroup prediction' finished... ---")


if __name__ == "__main__":
    main()
