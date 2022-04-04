import json
import math


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


def main():
    file = "test/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522_chrY/HG00096.mapped.ILLUMINA.bwa.GBR" \
           ".low_coverage.20120522_chrY.out"
    haplotype_dict = read_derived_haplotype_dict(file)
    read_backbone_groups()
    tree = Tree("Position_files/tree.json")
    get_most_likely_haplotype(tree, haplotype_dict)


def read_derived_haplotype_dict(file):
    # read all the derived lines and count the occurance of each branch
    haplotype_dict = {}
    with open(file) as f:
        f.readline()
        for line in f:
            _, _, _, branches, _, _, _, _, _, _, derived = line.strip().split("\t")
            for branch in branches.split(";"):
                if branch in haplotype_dict:
                    if derived == "D":
                        haplotype_dict[branch][0] += 1
                    haplotype_dict[branch][1] += 1
                else:
                    haplotype_dict[branch] = [0, 1]
                    if derived == "D":
                        haplotype_dict[branch][0] += 1
    return haplotype_dict


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


def get_most_likely_haplotype(tree, haplotype_dict):

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
                    qc3_score[0] += haplotype_dict[parent.name][0]
                    qc3_score[1] += haplotype_dict[parent.name][1]
            parent = parent.parent

        qc1_score = get_qc1_score(path, haplotype_dict)

        if haplotype_dict[node.name][1] == 0:
            qc2_score = 0
        else:
            qc2_score = haplotype_dict[node.name][0] / haplotype_dict[node.name][1]
        if qc3_score[1] == 0:
            qc3_score = 0
        else:
            qc3_score = qc3_score[0] / qc3_score[1]

        total_score = math.prod([qc1_score, qc2_score, qc3_score])

        if node.name in covered_node_scores and total_score <= covered_node_scores[node.name]:
            continue
        scores[node.name] = [qc1_score, qc2_score, qc3_score]
        tree_paths[node.name] = get_str_path(path)

        # make sure that less specific nodes with lower scores are not recorded
        for node_name in path:
            if node_name in covered_node_scores:
                if total_score > covered_node_scores[node_name]:
                    covered_node_scores[node_name] = total_score
            else:
                covered_node_scores[node_name] = total_score

    scores = sorted(scores.items(), key=lambda kv: (math.prod(kv[1]), tree.get(kv[0]).depth), reverse=True)
    for score in scores:
        total_score = math.prod(score[1])
        if total_score > 0:
            print(score[0], score[1], total_score)
            print(tree_paths[score[0]])


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
    for name, (derived, total) in intermediate_states.items():
        expected_state = expected_states[name]
        if expected_state == "A":
            score[0] += total - derived
        else:
            score[0] += derived
        score[1] += total
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


def write_outcome():
    pass


if __name__ == '__main__':
    main()
