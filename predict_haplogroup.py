#!/usr/bin/env python

"""
Script for Haplogroup prediction using the YFull tree

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import argparse
from pathlib import Path
from tree import Tree
from typing import Set, Dict, Iterator, List, Any, Union

BACKBONE_GROUPS: Set = set()
MAIN_HAPLO_GROUPS: Set = set()
QC1_SCORE_CACHE: Dict[str, float] = {}


class HgMarkersLinker:
    """Safe for a certain haplogroup if the number of ancestral and derived markers."""
    DERIVED: str = "D"
    ANCESTRAL: str = "A"
    UNDEFINED: str = "N"

    _ancestral_markers: Set[str]
    _derived_markers: Set[str]

    def __init__(self):
        self._ancestral_markers = set()
        self._derived_markers = set()

    def add(
        self,
        marker_name: str,
        state: str
    ):
        if state == self.DERIVED:
            self._derived_markers.add(marker_name)
        else:
            self._ancestral_markers.add(marker_name)

    def get_derived_markers(self) -> Set[str]:
        return self._derived_markers

    @property
    def nr_total(self) -> int:
        return len(self._ancestral_markers) + len(self._derived_markers)

    @property
    def nr_derived(self) -> int:
        return len(self._derived_markers)

    @property
    def nr_ancestral(self) -> int:
        return len(self._ancestral_markers)

    def get_state(self) -> str:
        """at least a fraction of 0.6 need to either derived or ancestral, otherwise the state can not be accurately
        determined and will be returned as undefined"""

        if self.nr_derived / self.nr_total >= 0.6:
            return self.DERIVED
        if self.nr_ancestral / self.nr_total >= 0.6:
            return self.ANCESTRAL
        return self.UNDEFINED


def main():
    """Main entry point for prediction script"""
    print("\tY-Haplogroup Prediction")
    namespace = get_arguments()
    in_folder = namespace.input
    output = namespace.outfile
    read_backbone_groups()
    final_table = []
    for folder in read_input_folder(in_folder):
        # make sure to reset this for each sample
        global QC1_SCORE_CACHE
        QC1_SCORE_CACHE = {}
        haplotype_dict = read_yleaf_out_file(folder / (folder.name + ".out"))
        tree = Tree("Hg_Prediction_tables/tree.json")
        best_haplotype_score = get_most_likely_haplotype(tree, haplotype_dict, namespace.minimum_score)
        add_to_final_table(final_table, haplotype_dict, best_haplotype_score, folder)
    write_final_table(final_table, output)
    print("--- Yleaf 'Y-Haplogroup prediction' finished... ---")


def get_arguments() -> argparse.Namespace:
    """Get the arguments provided by the user to this script"""
    parser = argparse.ArgumentParser(description="Erasmus MC: Genetic Identification\n Y-Haplogroup Prediction")

    parser.add_argument("-i", "--input", required=True,
                        help="Output file or path produced from Yleaf", metavar="FILE")

    parser.add_argument("-ms", "--minimum_score", help="Minimum score needed in order for a prediction to be considered"
                                                       "for inclusion in the final data (default=0.7).",
                        type=float, default=0.7)

    parser.add_argument("-o", "--outfile", required=True, help="Output file name", metavar="FILE")

    args = parser.parse_args()
    return args


def read_backbone_groups():
    """Read some basic data that is always needed"""
    global BACKBONE_GROUPS, MAIN_HAPLO_GROUPS
    with open("Hg_Prediction_tables/Intermediates.txt") as f:
        for line in f:
            if "~" in line:
                continue
            BACKBONE_GROUPS.add(line.strip())

    # add capitals A-Z to get all groups
    for nr in range(66, 85):
        MAIN_HAPLO_GROUPS.add(chr(nr))
    MAIN_HAPLO_GROUPS.add("A0-T")
    MAIN_HAPLO_GROUPS.add("A00")


def read_input_folder(
    folder_name: str
) -> Iterator[Path]:
    """Read all the folders present in the input folder. These folders are assumed to be created by Yleaf containing
     all the relevant information"""
    in_folder = Path(folder_name)
    for folder in in_folder.iterdir():
        if not folder.is_dir():
            continue
        yield folder


def read_yleaf_out_file(
    file: Union[Path, str]
) -> Dict[str, HgMarkersLinker]:
    """Read the full yleaf .out file and parse all lines into a dictionary keyed on haplogroups"""
    haplotype_dict = {}
    with open(file) as f:
        f.readline()
        for line in f:
            _, _, marker, haplogroup, _, _, _, _, _, _, state, _ = line.strip().split("\t")
            if haplogroup not in haplotype_dict:
                haplotype_dict[haplogroup] = HgMarkersLinker()
            haplotype_dict[haplogroup].add(marker, state)
    return haplotype_dict


def get_most_likely_haplotype(
    tree: Tree,
    haplotype_dict: Dict[str, HgMarkersLinker],
    treshold: float
) -> List[Any]:
    """Get the most specific haplogroup with a score above the treshold. The score is build up from 3 parts:
    part 1: This score indicates whether the predicted haplogroup follows the
        expected backbone of the haplogroup tree structure (i.e. if haplogroup E
        is predicted the markers defining: A0-T, A1, A1b, BT, CT, DE should be
        in the derived state, while other intermediate markers like: CF, F, GHIJK,
        etc, are expected to be in the ancestral state). The score is calculated by
        dividing the number of markers that show the expected state, by the sum
        of all intermediate markers. A score of 1 shows that all markers are in the
        expected state and indicates high confidence if the prediction of the correct
        broad haplogroup.
    part 2: This score indicates whether equivalent markers to the final haplogroup
        prediction were found in the ancestral state. I.e. if the final
        haplogroup is R1b1a1a2a2, there are two markers in the assay defining
        this haplogroup: Z2103 and Z2105, if both are found to be derived the
        part 2 score value will be 1. However if one is in the derived and the other in
        the ancestral state the QC2 value will be calculated as number of derived
        equivalent markers divided by the total number of equivalent markers, in
        this example the part 2 score value would be 0.5.
    part 3: This score indicates whether the predicted haplogroup follows the
        expected within-haplogroup tree structure. I.e. if the predicted haplogroup
        is O2a1c (O-JST002611), it is expected that markers defining:
        O2a1, O2a, O2 and O are also in the derived state. A score of 1 shows
        that all markers are in the expected state and indicates high confidence
        in the haplogroup prediction.

    First all possible haplogroups are sorted based on the depth in the tree. This is the order haplogroup scores will
    be calculated to determine if they are above the treshold.
    """
    sorted_depth_haplotypes = sorted(haplotype_dict.keys(), key=lambda k: tree.get(k).depth, reverse=True)
    covered_nodes = set()
    best_score = ["NA", "NA", 0, 0, 0, 0, 0]
    for haplotype_name in sorted_depth_haplotypes:
        node = tree.get(haplotype_name)

        # only record most specific nodes regardless of scores
        if node in covered_nodes:
            continue
        parent = node.parent
        path = [node.name]

        # caclulate score 3 already since this is most efficient
        qc3_score = [0, 0]  # matching, total

        # walk the tree back
        while parent is not None:
            path.append(parent.name)
            if parent.name in haplotype_dict:
                # in the same overal haplogroup
                if parent.name[0] == node.name[0] and parent.name != node.name[0]:
                    state = haplotype_dict[parent.name].get_state()
                    if state == HgMarkersLinker.DERIVED:
                        qc3_score[0] += 1

                    # in case it can not be decided what the state is ratio between 0.4 and  0.6
                    if state != HgMarkersLinker.UNDEFINED:
                        qc3_score[1] += 1
            parent = parent.parent

        qc1_score = get_qc1_score(path, haplotype_dict)
        if qc1_score < treshold:
            continue

        if haplotype_dict[node.name].nr_total == 0:
            qc2_score = 0
        else:
            qc2_score = haplotype_dict[node.name].nr_derived / haplotype_dict[node.name].nr_total
            if qc2_score < treshold:
                continue
        if qc3_score[1] == 0:
            qc3_score = 0
        else:
            qc3_score = qc3_score[0] / qc3_score[1]
            if qc3_score < treshold:
                continue

        total_score = product([qc1_score, qc2_score, qc3_score])

        # if above filter we found the hit
        if total_score > treshold:
            ancestral_children = get_ancestral_children(node, haplotype_dict)
            best_score = [node.name, ancestral_children, qc1_score, qc2_score, qc3_score, total_score, node.depth]
            break

        # make sure that less specific nodes are not recorded
        for node_name in path:
            covered_nodes.add(node_name)
    return best_score


def get_qc1_score(path, haplotype_dict):
    most_specific_backbone = None
    intermediate_states = {value: haplotype_dict[value] for value in BACKBONE_GROUPS if value in haplotype_dict}
    for value in path:
        if value in MAIN_HAPLO_GROUPS:
            most_specific_backbone = value
            break

    if most_specific_backbone is None:
        return 0

    global QC1_SCORE_CACHE
    # same backbone, same score
    if most_specific_backbone in QC1_SCORE_CACHE:
        return QC1_SCORE_CACHE[most_specific_backbone]
    else:
        expected_states = {}
        with open(f"Hg_Prediction_tables/{most_specific_backbone}_int.txt") as f:
            for line in f:
                if line == "":
                    continue
                name, state = line.strip().split("\t")
                expected_states[name] = {*state.split("/")}

    score = [0, 0]  # matching, total
    for name, marker_linker in intermediate_states.items():
        expected_possible_states = expected_states[name]
        state = marker_linker.get_state()
        if state in expected_possible_states:
            score[0] += 1
        if state != HgMarkersLinker.UNDEFINED:
            score[1] += 1
    if score[1] == 0:
        return 0

    # cache the score
    qc1_score = score[0] / score[1]
    QC1_SCORE_CACHE[most_specific_backbone] = qc1_score
    return qc1_score


def get_str_path(path):
    str_path = ""
    for value in path[:-1]:
        str_path = f"->{value}" + str_path
    str_path = "ROOT" + str_path
    return str_path


def get_ancestral_children(node, haplotype_dict):
    ancestral_children = []
    for name in node.children:
        if name in haplotype_dict:
            state = haplotype_dict[name].get_state()
            if state == HgMarkersLinker.ANCESTRAL:
                ancestral_children.append(name)
    return ancestral_children


def add_to_final_table(final_table, haplotype_dict, best_haplotype_scores, folder):
    total_reads, valid_markers = process_log_file(folder / (folder.name + ".log"))
    hg, ancestral_children, qc1, qc2, qc3, total, _ = best_haplotype_scores
    if hg == "NA":
        final_table.append([folder.name, hg, "", total_reads,
                            valid_markers, total, qc1, qc2, qc3])
        return
    marker_list = list(haplotype_dict[hg].get_derived_markers())
    # dont show all markers if too many
    if len(marker_list) > 2:
        marker_list = marker_list[:2] + ["etc."]
    if len(ancestral_children) > 0:
        ancestral_string = "x" + ','.join([f'{name.split("-")[1] if "-" in name else name}'
                                           for name in ancestral_children])
        final_table.append([folder.name, f"{hg}*({ancestral_string})", ';'.join(marker_list), total_reads,
                            valid_markers, total, qc1, qc2, qc3])
    else:
        final_table.append([folder.name, hg, ';'.join(marker_list), total_reads,
                            valid_markers, total, qc1, qc2, qc3])


def process_log_file(log_file):
    total_reads = None
    valid_markers = None

    with open(log_file) as f:
        for line in f:
            if line.startswith("Total of reads"):
                total_reads = int(line.replace("Total of reads:", "").strip())
            elif line.startswith("Markers with haplogroup information"):
                valid_markers = int(line.replace("Markers with haplogroup information:", "").strip())
    return total_reads, valid_markers


def write_final_table(final_table, out_file):
    # sort for each sample based on QC-score
    final_table.sort(key=lambda values: (values[0], values[5]), reverse=True)

    header = "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n"
    with open(out_file, "w") as f:
        f.write(header)
        for line in final_table:
            f.write('\t'.join(map(str, line)) + "\n")


def product(values):
    total = 1
    for value in values:
        total *= value
    return total


if __name__ == '__main__':
    main()
