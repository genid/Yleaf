import json
import math
import argparse
from pathlib import Path

BACKBONE_GROUPS = set()
MAIN_HAPLO_GROUPS = set()
QC1_SCORE_CACHE = {}


class Tree:
    __ROOT_KEY = 'ROOT (Y-Chromosome "Adam")'

    def __init__(self, file):
        self.node_mapping = {}
        self._construct_tree(file)

    def _construct_tree(self, file):
        with open(file) as f:
            tree_dict = json.load(f)

        # base has no parent
        base_node = Node(self.__ROOT_KEY, None, 0)
        self.node_mapping[self.__ROOT_KEY] = base_node
        self._recursive_read(tree_dict, base_node, 0)

    def _recursive_read(self, tree_dict, node, depth):
        depth += 1
        for node_name in tree_dict[node.name]:
            new_node = Node(node_name, node, depth)
            self.node_mapping[node_name] = new_node
            self._recursive_read(tree_dict, new_node, depth)

    def get(self, node_name):
        return self.node_mapping[node_name]


class Node:
    def __init__(self, name, parent_node, depth):
        self.name = name
        self.parent = parent_node
        self.depth = depth


class HgMarkersLinker:
    DERIVED = "D"
    ANCESTRAL = "A"

    def __init__(self):
        self._ancestral_markers = set()
        self._derived_markers = set()

    def add(self, marker_name, state):
        if state == self.DERIVED:
            self._derived_markers.add(marker_name)
        else:
            self._ancestral_markers.add(marker_name)

    def get_derived_markers(self):
        return self._derived_markers

    @property
    def nr_total(self):
        return len(self._ancestral_markers) + len(self._derived_markers)

    @property
    def nr_derived(self):
        return len(self._derived_markers)

    @property
    def nr_ancestral(self):
        return len(self._ancestral_markers)


def main():
    print("\tY-Haplogroup Prediction")
    namespace = get_arguments()
    in_folder = namespace.input
    output = namespace.outfile
    read_backbone_groups()
    final_table = []
    for folder in read_input_folder(in_folder):
        haplotype_dict = read_derived_haplotype_dict(folder / (folder.name + ".out"))
        tree = Tree("Position_files/tree.json")
        best_haplotype_scores = get_most_likely_haplotype(tree, haplotype_dict)
        add_to_final_table(final_table, haplotype_dict, best_haplotype_scores, folder)
    write_final_table(final_table, output)
    print("--- Yleaf 'Y-Haplogroup prediction' finished... ---")


def get_arguments():
    parser = argparse.ArgumentParser(description="Erasmus MC: Genetic Identification\n Y-Haplogroup Prediction")

    parser.add_argument("-i", "--input", required=True,
                        help="Output file or path produced from Yleaf", metavar="FILE")

    parser.add_argument("-o", "--outfile", required=True, help="Output file name", metavar="FILE")

    args = parser.parse_args()
    return args


def read_backbone_groups():
    global BACKBONE_GROUPS, MAIN_HAPLO_GROUPS
    with open("Hg_Prediction_tables/Intermediates.txt") as f:
        for line in f:
            if "~" in line:
                continue
            BACKBONE_GROUPS.add(line.strip())

    # add capitals A-Z to get all groups
    for nr in range(65, 85):
        MAIN_HAPLO_GROUPS.add(chr(nr))


def read_input_folder(folder_name):
    in_folder = Path(folder_name)
    for folder in in_folder.iterdir():
        if not folder.is_dir():
            continue
        yield folder


def read_derived_haplotype_dict(file):
    # read all the derived lines and count the occurance of each branch
    haplotype_dict = {}
    with open(file) as f:
        f.readline()
        for line in f:
            _, _, marker, branch, _, _, _, _, _, _, state = line.strip().split("\t")
            if branch not in haplotype_dict:
                haplotype_dict[branch] = HgMarkersLinker()
            haplotype_dict[branch].add(marker, state)
    return haplotype_dict


def get_most_likely_haplotype(tree, haplotype_dict, treshold=0.7):

    sorted_depth_haplotypes = sorted(haplotype_dict.keys(), key=lambda k: tree.get(k).depth, reverse=True)
    covered_node_scores = {}
    scores = {}
    tree_paths = {}
    for haplotype_name in sorted_depth_haplotypes:
        node = tree.get(haplotype_name)
        parent = node.parent
        path = [node.name]

        # caclulate score 3 already since this is most efficient
        qc3_score = [0, 0]  # matching, total

        # walk the tree back
        while parent is not None:
            path.append(parent.name)
            # at least one of the snps is in derived state and in the same overal haplogroup
            if parent.name in haplotype_dict:
                if parent.name[0] == node.name[0] and parent.name != node.name[0]:
                    qc3_score[0] += haplotype_dict[parent.name].nr_derived
                    qc3_score[1] += haplotype_dict[parent.name].nr_total
            parent = parent.parent

        qc1_score = get_qc1_score(path, haplotype_dict)

        if haplotype_dict[node.name].nr_total == 0:
            qc2_score = 0
        else:
            qc2_score = haplotype_dict[node.name].nr_derived / haplotype_dict[node.name].nr_total
        if qc3_score[1] == 0:
            qc3_score = 0
        else:
            qc3_score = qc3_score[0] / qc3_score[1]

        total_score = math.prod([qc1_score, qc2_score, qc3_score])

        if node.name in covered_node_scores and total_score <= covered_node_scores[node.name]:
            continue
        # only record these scores
        if total_score > treshold:
            scores[node.name] = [qc1_score, qc2_score, qc3_score]
            tree_paths[node.name] = get_str_path(path)

        # make sure that less specific nodes with lower scores are not recorded
        for node_name in path:
            if node_name in covered_node_scores:
                if total_score > covered_node_scores[node_name]:
                    covered_node_scores[node_name] = total_score
            else:
                covered_node_scores[node_name] = total_score

    return scores


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
                expected_states[name] = state

    score = [0, 0]  # matching, total
    for name, marker_linker in intermediate_states.items():
        expected_state = expected_states[name]
        if expected_state == "A":
            score[0] += marker_linker.nr_ancestral
        else:
            score[0] += marker_linker.nr_derived
        score[1] += marker_linker.nr_total
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


def add_to_final_table(final_table, haplotype_dict, best_haplotype_scores, folder):
    total_reads, valid_markers = process_log_file(folder / (folder.name + ".log"))
    for hg, score_values in best_haplotype_scores.items():
        final_table.append([folder.name, hg, ';'.join(haplotype_dict[hg].get_derived_markers()), total_reads,
                            valid_markers, math.prod(score_values), *score_values])


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


if __name__ == '__main__':
    main()
