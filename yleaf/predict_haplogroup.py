#!/usr/bin/env python

"""
Script for Haplogroup prediction using the YFull tree

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import argparse
import json
import multiprocessing
from functools import partial
from pathlib import Path
from typing import Set, Dict, Iterator, List, Any, Union, Tuple
import logging

from yleaf import yleaf_constants
from yleaf.tree import Tree, Node

BACKBONE_GROUPS: Set = set()
MAIN_HAPLO_GROUPS: Set = set()
QC1_SCORE_CACHE: Dict[str, float] = {}
ACTIVE_TREE: str = "yfull"  # tracks which tree is active for QC routing
_WORKER_TREE: 'Tree' = None  # cached per worker, set by _init_predict_worker

DEFAULT_MIN_SCORE = 0.95
JSON_REPORT_SCHEMA_VERSION = 1

LOG = logging.getLogger("yleaf_logger")


def _init_predict_worker(tree_file: Path):
    global _WORKER_TREE
    _WORKER_TREE = Tree(tree_file)


class HgMarkersLinker:
    """Safe for a certain haplogroup if the number of ancestral and derived markers."""
    DERIVED: str = "D"
    ANCESTRAL: str = "A"
    UNDEFINED: str = "N"

    _ancestral_markers: Set[str]
    _derived_markers: Set[str]
    _d: int
    _a: int

    def __init__(self):
        self._ancestral_markers = set()
        self._derived_markers = set()
        self._d = 0
        self._a = 0

    def add(
        self,
        marker_name: str,
        state: str
    ):
        if state == self.DERIVED:
            self._derived_markers.add(marker_name)
            self._d += 1
        else:
            self._ancestral_markers.add(marker_name)
            self._a += 1

    def get_derived_markers(self) -> Set[str]:
        return self._derived_markers

    @property
    def nr_total(self) -> int:
        return self._d + self._a

    @property
    def nr_derived(self) -> int:
        return self._d

    @property
    def nr_ancestral(self) -> int:
        return self._a

    def get_state(self) -> str:
        """at least a fraction of 0.6 need to either derived or ancestral, otherwise the state can not be accurately
        determined and will be returned as undefined"""
        d, a = self._d, self._a
        total = d + a
        if total == 0:
            return self.UNDEFINED
        if d / total >= 0.6:
            return self.DERIVED
        if a / total >= 0.6:
            return self.ANCESTRAL
        return self.UNDEFINED


def main_predict_haplogroup(
    namespace: argparse.Namespace,
    backbone_groups: Set,
    main_haplo_groups: Set,
    folder: Path,
    active_tree: str = "yfull",
    out_file_suffix: str = ".out",
):
    # Re-assign globals so spawned worker processes (which don't inherit parent
    # state on platforms using the "spawn" start method, e.g. macOS) have the
    # required data available (fixes issue #34).
    global BACKBONE_GROUPS, MAIN_HAPLO_GROUPS, QC1_SCORE_CACHE, ACTIVE_TREE
    BACKBONE_GROUPS = backbone_groups
    MAIN_HAPLO_GROUPS = main_haplo_groups
    ACTIVE_TREE = active_tree
    # make sure to reset this for each sample
    QC1_SCORE_CACHE = {}

    if folder.name == "filtered_vcf_files":
        return [None, None, None]

    out_path = folder / (folder.name + out_file_suffix)
    try:
        haplotype_dict = read_yleaf_out_file(out_path)
    except FileNotFoundError:
        LOG.warning(f"WARNING: failed to find .out file from yleaf run for sample {folder.name}. This sample will"
                    " be skipped.")
        return [None, None, None]
    LOG.debug(f"[diag] {folder.name}: read {len(haplotype_dict)} haplogroups from {out_path.name}; "
              f"ACTIVE_TREE={ACTIVE_TREE} BACKBONE={len(BACKBONE_GROUPS)} MAIN_HG={len(MAIN_HAPLO_GROUPS)}")
    if _WORKER_TREE is not None:
        tree = _WORKER_TREE
    else:
        tree_file = getattr(namespace, 'tree_file', None) or (yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.TREE_FILE)
        tree = Tree(tree_file)
    best_haplotype_score = get_most_likely_haplotype(tree, haplotype_dict, namespace.minimum_score)
    LOG.debug(f"[diag] {folder.name}: best={best_haplotype_score[0]}")
    return [haplotype_dict, best_haplotype_score, folder]


def main(namespace: argparse.Namespace = None):
    """Main entry point for prediction script"""
    LOG.info("Starting haplogroup prediction...")
    if namespace is None:
        namespace = get_arguments()
    in_folder = namespace.input
    output = namespace.outfile
    threads = namespace.threads
    tree_file = getattr(namespace, 'tree_file', None)
    out_file_suffix = getattr(namespace, 'out_file_suffix', '.out')
    info_file_suffix = getattr(namespace, 'info_file_suffix', '.info')
    report_json = getattr(namespace, 'report_json', None)
    read_backbone_groups(tree_file)
    final_table = []
    final_records = [] if report_json else None

    folders = list(read_input_folder(in_folder))
    total = len(folders)
    LOG.info(f"[PROGRESS] prediction 0/{total}")
    effective_tree_file = tree_file or (yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.TREE_FILE)
    with multiprocessing.Pool(processes=threads,
                               initializer=_init_predict_worker,
                               initargs=(effective_tree_file,)) as p:
        for i, result in enumerate(
            p.imap_unordered(partial(main_predict_haplogroup, namespace, BACKBONE_GROUPS, MAIN_HAPLO_GROUPS,
                                     active_tree=ACTIVE_TREE, out_file_suffix=out_file_suffix),
                             folders),
            1,
        ):
            haplotype_dict, best_haplotype_score, folder = result
            if haplotype_dict is not None:
                add_to_final_table(final_table, haplotype_dict, best_haplotype_score, folder, info_file_suffix)
                if final_records is not None:
                    add_to_json_records(final_records, haplotype_dict, best_haplotype_score, folder, info_file_suffix)
            LOG.info(f"[PROGRESS] prediction {i}/{total}")

    write_final_table(final_table, output)
    if final_records is not None:
        write_json_report(final_records, report_json, tree_name=ACTIVE_TREE)
    LOG.debug("Finished haplogroup prediction")


def get_arguments() -> argparse.Namespace:
    """Get the arguments provided by the user to this script"""
    parser = argparse.ArgumentParser(description="Erasmus MC: Genetic Identification\n Y-Haplogroup Prediction")

    parser.add_argument("-i", "--input", required=True,
                        help="Output file or path produced from Yleaf", metavar="FILE")

    parser.add_argument("-ms", "--minimum_score", help="Minimum score needed in order for a prediction to be considered"
                                                       "for inclusion in the final data (default=0.95).",
                        type=float, default=DEFAULT_MIN_SCORE)

    parser.add_argument("-o", "--outfile", required=True, help="Output file name", metavar="FILE")

    parser.add_argument("-t", "--threads", help="Number of threads to use (default=1).", type=int, default=1)

    args = parser.parse_args()
    return args


def read_backbone_groups(tree_file: Path = None):
    """Read some basic data that is always needed"""
    global BACKBONE_GROUPS, MAIN_HAPLO_GROUPS, ACTIVE_TREE
    # Reset per-tree globals so sequential multi-tree calls don't accumulate stale entries
    BACKBONE_GROUPS = set()
    MAIN_HAPLO_GROUPS = set()
    hg_folder = yleaf_constants.HG_PREDICTION_FOLDER

    is_ftdna = tree_file is not None and str(tree_file).endswith(yleaf_constants.FTDNA_TREE_FILE)
    is_isogg = tree_file is not None and str(tree_file).endswith(yleaf_constants.ISOGG_TREE_FILE)
    is_yfull_v10 = tree_file is not None and str(tree_file).endswith(yleaf_constants.YFULL_V10_TREE_FILE)

    if is_ftdna:
        ACTIVE_TREE = yleaf_constants.TREE_FTDNA
        major_tables_dir = hg_folder / "ftdna_major_tables"
    elif is_isogg:
        ACTIVE_TREE = yleaf_constants.TREE_ISOGG
        major_tables_dir = hg_folder / "isogg_major_tables"
    elif is_yfull_v10:
        ACTIVE_TREE = yleaf_constants.TREE_YFULL_V10
        major_tables_dir = hg_folder / "yfull_v10_major_tables"
    else:
        ACTIVE_TREE = yleaf_constants.TREE_YFULL
        major_tables_dir = hg_folder / "major_tables"

    with open(major_tables_dir / "Intermediates.txt") as f:
        for line in f:
            if "~" in line:
                continue
            name = line.strip()
            if name:
                BACKBONE_GROUPS.add(name)
                MAIN_HAPLO_GROUPS.add(name)

    if ACTIVE_TREE in (yleaf_constants.TREE_YFULL, yleaf_constants.TREE_YFULL_V10):
        # YFull (v10 and v14) also needs these intermediate nodes in MAIN_HAPLO_GROUPS
        major_tree_list = ['A00', 'A00a', 'A00b', 'A00c', 'A0-T', 'A1', 'A1b', 'A1b1', 'BT', 'CT', 'CF', 'F', 'F4', 'F2', 'F3', 'GHIJK', 'G', 'HIJK', 'H', 'IJK', 'IJ', 'I', 'I2', 'I1', 'J', 'J1', 'J2', 'K', 'K2', 'K2d', 'K2c', 'K2b', 'P', 'R', 'R2', 'R1', 'R1b', 'R1a', 'Q', 'K2b1', 'S', 'M', 'NO', 'O', 'N', 'LT', 'L', 'T', 'F1', 'C', 'DE', 'D', 'E', 'B', 'A1a', 'A0']
        for major_hg in major_tree_list:
            MAIN_HAPLO_GROUPS.add(major_hg)


def read_input_folder(
    folder_name: str
) -> Iterator[Path]:
    """Read all the folders present in the input folder. These folders are assumed to be created by Yleaf containing
     all the relevant information"""
    in_folder = Path(folder_name)
    for folder in in_folder.iterdir():
        if not folder.is_dir():
            continue
        if folder.name == yleaf_constants.FASTQ_BAM_FILE_FOLDER:
            # generated for efficiency during fastq, needs to be ignored
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
) -> Tuple[str, str, int, int, int, int, int]:
    """Get the most specific haplogroup with a score above the treshold. The score is build up from 3 parts:
    QC1: This score indicates whether the predicted haplogroup follows the
        expected backbone of the haplogroup tree structure (i.e. if haplogroup E
        is predicted the markers defining: A0-T, A1, A1b, BT, CT, DE should be
        in the derived state, while other intermediate markers like: CF, F, GHIJK,
        etc, are expected to be in the ancestral state). The score is calculated by
        dividing the number of markers that show the expected state, by the sum
        of all intermediate markers. A score of 1 shows that all markers are in the
        expected state and indicates high confidence if the prediction of the correct
        broad haplogroup.
    QC2: This score indicates whether equivalent markers to the final haplogroup
        prediction were found in the ancestral state. I.e. if the final
        haplogroup is R1b1a1a2a2, there are two markers in the assay defining
        this haplogroup: Z2103 and Z2105, if both are found to be derived the
        part 2 score value will be 1. However if one is in the derived and the other in
        the ancestral state the QC2 value will be calculated as number of derived
        equivalent markers divided by the total number of equivalent markers, in
        this example the part 2 score value would be 0.5.
    QC3: This score indicates whether the predicted haplogroup follows the
        expected within-haplogroup tree structure. I.e. if the predicted haplogroup
        is O2a1c (O-JST002611), it is expected that markers defining:
        O2a1, O2a, O2 and O are also in the derived state. A score of 1 shows
        that all markers are in the expected state and indicates high confidence
        in the haplogroup prediction.

    First all possible haplogroups are sorted based on the depth in the tree. This is the order haplogroup scores will
    be calculated to determine if they are above the treshold. As soon as the treshold is reached the function returns.
    """
    sorted_depth_haplotypes = sorted(haplotype_dict.keys(), key=lambda k: tree.get(k).depth, reverse=True)
    covered_nodes = set()
    best_score = ("NA", "NA", "NA", "NA", "NA", "NA", "NA")

    for haplotype_name in sorted_depth_haplotypes:
        node = tree.get(haplotype_name)

        # only record most specific nodes regardless of scores
        if node in covered_nodes:
            continue

        # QC2 pre-check: skip path building for haplogroups with no/insufficient derived markers.
        # Profiling shows 99.9 % of haplogroups fail here, so this avoids the expensive tree walk
        # for virtually all candidates.
        h = haplotype_dict[haplotype_name]
        _d, _a = h._d, h._a
        _total = _d + _a
        if _total == 0:
            qc2_score = 0
        else:
            qc2_score = _d / _total
        if qc2_score < treshold:
            continue

        node_name = node.name
        node_first_char = node_name[0]
        parent = node.parent
        path = [node_name]

        # caclulate score 3 already since this is most efficient
        qc3_score = [0, 0]  # matching, total

        # walk the tree back
        while parent is not None:
            pname = parent.name
            path.append(pname)
            if pname in haplotype_dict:
                # in the same overal haplogroup
                if node_first_char in pname and pname != node_first_char:
                    state = haplotype_dict[pname].get_state()
                    if state == HgMarkersLinker.DERIVED:
                        qc3_score[0] += 1

                    # in case it can not be decided what the state is ratio between 0.4 and 0.6
                    if state != HgMarkersLinker.UNDEFINED:
                        qc3_score[1] += 1

            parent = parent.parent

        qc1_score = get_qc1_score(path, haplotype_dict)

        # if any of the scores are below treshold, the total can not be above so ignore
        if qc1_score < treshold:
            continue

        # QC2 already computed above; only repeat the threshold check
        if qc2_score < treshold:
            continue
        if qc3_score[1] == 0:
            # No intermediate ancestor data — neutral (no contradicting evidence).
            # FTDNA keeps 0 to block false positives at low coverage (see CLAUDE.md).
            qc3_score = 0 if ACTIVE_TREE == yleaf_constants.TREE_FTDNA else 1.0
        else:
            qc3_score = qc3_score[0] / qc3_score[1]
            if qc3_score < treshold:
                continue

        total_score = product([qc1_score, qc2_score, qc3_score])

        # if above filter we found the hit
        if total_score >= treshold:
            ancestral_children = get_ancestral_children(node, haplotype_dict, tree)
            best_score = (node.name, ancestral_children, qc1_score, qc2_score, qc3_score, total_score, node.depth)
            break
            
        # else, no hit is found, but still report Quality Scores
        else:
            best_score = ("NA", "NA", qc1_score, qc2_score, qc3_score, total_score, "NA")
            
        # make sure that less specific nodes are not recorded
        for node_name in path:
            covered_nodes.add(node_name)
    return best_score


def get_qc1_score(
    path: List[str],
    haplotype_dict: Dict[str, HgMarkersLinker]
) -> float:
    """Get the first quality score as described in the get_most_likely_haplotype function"""
    global QC1_SCORE_CACHE
    most_specific_backbone = None
    hg_folder = yleaf_constants.HG_PREDICTION_FOLDER

    intermediate_states = {value: haplotype_dict[value] for value in BACKBONE_GROUPS if value in haplotype_dict}
    for value in path:
        if value in MAIN_HAPLO_GROUPS:
            most_specific_backbone = value
            break

    if most_specific_backbone is None:
        return 0

    if most_specific_backbone in QC1_SCORE_CACHE:
        return QC1_SCORE_CACHE[most_specific_backbone]

    if ACTIVE_TREE == yleaf_constants.TREE_FTDNA:
        tables_subdir = "ftdna_major_tables"
    elif ACTIVE_TREE == yleaf_constants.TREE_ISOGG:
        tables_subdir = "isogg_major_tables"
    else:
        tables_subdir = "major_tables"

    expected_states = {}
    with open(f"{hg_folder}/{tables_subdir}/{most_specific_backbone}_int.txt") as f:
            for line in f:
                if line == "":
                    continue
                name, state = line.strip().split("\t")
                expected_states[name] = {*state.split("/")}

    score = [0, 0]  # matching, total
    if ACTIVE_TREE == yleaf_constants.TREE_FTDNA:
        # All 20 backbone groups are always in the denominator.
        # Absent or UNDEFINED groups are neutral (no penalty); only clear mismatches reduce the score.
        for name, expected_possible_states in expected_states.items():
            score[1] += 1
            if name not in intermediate_states:
                score[0] += 1  # no data — neutral
            else:
                state = intermediate_states[name].get_state()
                if state == HgMarkersLinker.UNDEFINED or state in expected_possible_states:
                    score[0] += 1  # matches or undetermined — no penalty
    else:
        for name, marker_linker in intermediate_states.items():
            expected_possible_states = expected_states[name]
            state = marker_linker.get_state()
            if state in expected_possible_states:
                score[0] += 1
            if state != HgMarkersLinker.UNDEFINED:
                score[1] += 1
    if score[1] == 0:
        return 1.0  # no backbone observed — neutral, let QC2/QC3 decide

    # cache the score
    qc1_score = score[0] / score[1]
    QC1_SCORE_CACHE[most_specific_backbone] = qc1_score
    return qc1_score


def get_ancestral_children(
    node: Node,
    haplotype_dict: Dict[str, HgMarkersLinker],
    tree: Tree
) -> List[str]:
    """Get the names of the first nodes that are ancestral children of the given node do this until no nodes are left
    or all ancestral children are found"""
    ancestral_children = []
    # make sure to not modify the original
    to_cover = [node]
    while len(to_cover) > 0:
        node = to_cover.pop()
        for name in node.children:
            if name in haplotype_dict:
                state = haplotype_dict[name].get_state()
                if state == HgMarkersLinker.ANCESTRAL:
                    ancestral_children.append(name)
                    continue
            to_cover.append(tree.get(name))
    return ancestral_children


def add_to_final_table(
    final_table: List[List[Any]],
    haplotype_dict: Dict[str, HgMarkersLinker],
    best_haplotype_scores: Tuple[str, str, int, int, int, int, int],
    folder: Path,
    info_file_suffix: str = ".info"
):
    """Add the score to the final table and some additional information from the log file generated by Yleaf"""
    total_reads, valid_markers = process_info_file(folder / (folder.name + info_file_suffix))
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
        ancestral_string = "x" + ','.join(ancestral_children)
        final_table.append([folder.name, f"{hg}*({ancestral_string})", ';'.join(marker_list), total_reads,
                            valid_markers, total, qc1, qc2, qc3])
    else:
        final_table.append([folder.name, hg, ';'.join(marker_list), total_reads,
                            valid_markers, total, qc1, qc2, qc3])


def process_info_file(
    info_file: Path
) -> Tuple[Union[int, str], Union[int, str]]:
    """Read the log file produced by Yleaf for some additional information in the output"""
    total_reads = "NA"
    valid_markers = "NA"

    try:
        with open(info_file) as f:
            for line in f:
                if line.startswith("Total of mapped reads:"):
                    total_reads = line.replace("Total of mapped reads:", "").strip()
                elif line.startswith("Markers with haplogroup information"):
                    valid_markers = line.replace("Markers with haplogroup information:", "").strip()
    except FileNotFoundError:
        LOG.warning("WARNING: failed to find .info file from yleaf run. This information is not critical but there"
                    " will be some missing values in the output.")
    return total_reads, valid_markers


def write_final_table(
    final_table: List[List[Any]],
    out_file: Union[Path, str]
):
    # sort for each sample based on QC-score
    final_table.sort(key=lambda values: (values[0], values[5]), reverse=True)

    header = "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n"
    with open(out_file, "w") as f:
        f.write(header)
        for line in final_table:
            f.write('\t'.join(map(str, line)) + "\n")


def _maybe_int(val: Any) -> Union[int, None]:
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def _maybe_float(val: Any) -> Union[float, None]:
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def add_to_json_records(
    records: List[Dict],
    haplotype_dict: Dict[str, HgMarkersLinker],
    best_haplotype_scores: Tuple,
    folder: Path,
    info_file_suffix: str = ".info"
):
    """Add a structured per-sample record for the JSON report."""
    total_reads, valid_markers = process_info_file(folder / (folder.name + info_file_suffix))
    hg, ancestral_children, qc1, qc2, qc3, total, _ = best_haplotype_scores

    if hg == "NA":
        full_hg = "NA"
        hg_base = "NA"
        markers = []
    else:
        markers = sorted(haplotype_dict[hg].get_derived_markers()) if hg in haplotype_dict else []
        if ancestral_children:
            ancestral_string = "x" + ','.join(ancestral_children)
            full_hg = f"{hg}*({ancestral_string})"
        else:
            full_hg = hg
        hg_base = hg

    records.append({
        "sample_id": folder.name,
        "haplogroup": full_hg,
        "haplogroup_base": hg_base,
        "excluded_subclades": list(ancestral_children) if hg != "NA" else [],
        "markers": markers,
        "total_reads": _maybe_int(total_reads),
        "valid_markers": _maybe_int(valid_markers),
        "qc_score": _maybe_float(total),
        "qc1": _maybe_float(qc1),
        "qc2": _maybe_float(qc2),
        "qc3": _maybe_float(qc3),
    })


def write_json_report(
    records: List[Dict],
    out_file: Union[Path, str],
    tree_name: str = "yfull"
):
    """Write a schema-versioned JSON report alongside the TSV."""
    from yleaf import __version__
    report = {
        "schema_version": JSON_REPORT_SCHEMA_VERSION,
        "yleaf_version": __version__,
        "tree": tree_name,
        "samples": sorted(records, key=lambda r: r["sample_id"]),
    }
    with open(out_file, "w") as f:
        json.dump(report, f, indent=2)


def product(
    values: List[float]
) -> float:
    """Take the product of a list of floats. The math.prod function is python 3.8 or above"""
    total = 1
    for value in values:
        total *= value
    return total


if __name__ == '__main__':
    main()
