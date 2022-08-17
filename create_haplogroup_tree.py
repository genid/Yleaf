from typing import List, Dict, Tuple
import graphviz
from collections import defaultdict
import argparse

from tree import Tree


def main():
    namespace = get_arguments()
    haplogroups, sample_mapping = read_input_file(namespace.input)
    partial_haplogroup_dict = haplogroup_tree_dict(haplogroups)
    make_dendrogram(partial_haplogroup_dict, sample_mapping, namespace.outfile)


def get_arguments() -> argparse.Namespace:
    """Get the arguments provided by the user to this script"""
    parser = argparse.ArgumentParser(description="Erasmus MC: Genetic Identification\n Haplogroup tree creation")

    parser.add_argument("-i", "--input", required=True,
                        help="prediction file (.hg file) containing the predicted haplogroups for different samples",
                        metavar="FILE")

    parser.add_argument("-o", "--outfile", required=True, help="Output file name", metavar="FILE")

    args = parser.parse_args()
    return args


def read_input_file(
    file: str
) -> Tuple[List[str], Dict[str, List[str]]]:
    haplogroups = []
    sample_mapping = defaultdict(list)
    with open(file) as f:
        f.readline()
        for line in f:
            values = line.split("\t")
            sample, haplogroup = values[:2]
            if haplogroup == 'NA':
                continue
            haplogroup = haplogroup.split("*")[0]
            haplogroups.append(haplogroup)
            sample_mapping[haplogroup].append(sample)
    return haplogroups, sample_mapping


def haplogroup_tree_dict(
    haplogroups: List[str]
) -> Dict[str, List[str]]:
    tree = Tree("Hg_Prediction_tables/tree.json")
    partial_tree_dict = defaultdict(list)
    for name in haplogroups:
        path = []
        node = tree.get(name)
        while node is not None:
            path.append(node.name)
            node = node.parent
        parent = path.pop()
        while len(path) > 0:
            child = path.pop()
            partial_tree_dict[parent].append(child)
            parent = child
    return partial_tree_dict


def make_dendrogram(
    partial_tree_dict: Dict[str, List[str]],
    sample_mapping: Dict[str, List[str]],
    output_file: str
):
    dot = graphviz.Digraph()
    dot.attr(ratio="compress")
    covered_nodes = set()
    covered_edges = set()
    for parent, children in partial_tree_dict.items():
        if parent not in covered_nodes:
            dot.node(parent, parent, shape="box")
            covered_nodes.add(parent)
        for child in children:
            if child not in covered_nodes:
                dot.node(child, child, shape="box")
                covered_nodes.add(child)
            if (parent, child) not in covered_edges:
                dot.edge(parent, child, weight="1")
                covered_edges.add((parent, child))
                if child in sample_mapping:
                    dot.node('\n'.join(sample_mapping[child]), '\n'.join(sample_mapping[child]), shape="box",
                             style='filled', fillcolor="green")
                    dot.edge(child, '\n'.join(sample_mapping[child]), weight="10")
    dot.render(output_file, view=False)


if __name__ == '__main__':
    main()
