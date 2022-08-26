#!/usr/bin/env python

"""
Script for transforming the ISOGG csv table into tables that can be used by Yleaf

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import argparse
from typing import List, Tuple


def main():
    namespace = get_arguments()
    hg19_data, hg38_data = read_isogg_data(namespace.input)
    write_table(hg19_data, hg38_data, namespace.out_name)


def get_arguments() -> argparse.Namespace:
    """Get the arguments provided by the user to this script"""
    parser = argparse.ArgumentParser(description="Transform the data table that can be downloaded from the ISOGG "
                                                 "website into tables that can be read by Yleaf")

    parser.add_argument("-i", "--input", required=True,
                        help="ISOGG data table, this can be downloaded from the https://isogg.org/ website",
                        metavar="FILE")

    parser.add_argument("-o", "--out_name", required=True, help="First part of the output file name. Always 2 files"
                                                                " will be produced. One for hg19 and one for hg38",
                        metavar="FILE")

    args = parser.parse_args()
    return args


def read_isogg_data(
    infile: str
) -> Tuple[List[Tuple], List[Tuple]]:
    """Convert the ISOGG table into two tables for 2 versions of the human genome. """
    skipped_isog = 0
    hg19_data = []
    hg38_data = []
    with open(infile) as f:
        for line in f:
            try:
                name, branch, alternate_names, _, hg19_location, hg38_location, ancestor_derived = \
                    line.strip().replace(" ", "").split(",")
            except ValueError:
                skipped_isog += 1
                continue
            try:
                ancestor, derived = ancestor_derived.split("->")
            except ValueError:
                # invalid change
                skipped_isog += 1
                continue
            if ".." in hg19_location:
                hg19_location = hg19_location.split("..")[0]
            if ".." in hg38_location:
                hg38_location = hg38_location.split("..")[0]

            # filter all these since they are likely faulty
            if "(" in branch or '[' in branch or branch == "#REF!" or branch == '':
                continue

            alternate_names = alternate_names.replace('"', "").split(";")
            if name != "":
                if hg19_location != "None":
                    hg19_data.append(("chry", name, branch, hg19_location, f"{ancestor}->{derived}", ancestor, derived))
                if hg38_location != "None":
                    hg38_data.append(("chry", name, branch, hg38_location, f"{ancestor}->{derived}", ancestor, derived))

            for aname in alternate_names:
                if aname == "":
                    continue
                if hg19_location != "None":
                    hg19_data.append(("chry", aname, branch, hg19_location, f"{ancestor}->{derived}", ancestor,
                                      derived))
                if hg38_location != "None":
                    hg38_data.append(("chry", aname, branch, hg38_location, f"{ancestor}->{derived}", ancestor,
                                      derived))
    print(f"Skipped {skipped_isog} isogg rows because format is inconsistent")
    hg19_data.sort(key=lambda values: values[2])
    hg38_data.sort(key=lambda values: values[2])
    return hg19_data, hg38_data


def write_table(
    hg19_data: List[Tuple],
    hg38_data: List[Tuple],
    out_name: str
):
    """Write the tables into tsv format but keep the .txt to keep consistent with old versions"""
    with open(f"{out_name}_hg19.txt", "w") as f:
        for line in hg19_data:
            f.write('\t'.join(line) + '\n')

    with open(f"{out_name}_hg38.txt", "w") as f:
        for line in hg38_data:
            f.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()
