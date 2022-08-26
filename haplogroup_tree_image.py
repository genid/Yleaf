#!/usr/bin/env python

"""
Script for drawing haplogroups in the haplogroup tree. Is designed to be used with Yleaf but can also be used with
Yfull haplogroups

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import os
from typing import List, Dict, Tuple, Set
import graphviz
from collections import defaultdict
import argparse
import logging

from tree import Tree

LOG = logging.getLogger("yleaf_logger")


def main(namespace: argparse.Namespace = None):
    if namespace is None:
        namespace = get_arguments()
    LOG.info("Starting with drawing haplogroups")
    haplogroups, sample_mapping = read_input_file(namespace.input)
    if len(haplogroups) == 0:
        LOG.warning("No haplogroups found in provided input file.")
        return
    add_main_haplogroups(haplogroups)
    partial_haplogroup_dict = haplogroup_tree_dict(haplogroups)
    if namespace.collapse_mode:
        edge_mapping = collapse_tree_dict(partial_haplogroup_dict, sample_mapping)
    else:
        edge_mapping = {}
    make_dendrogram(partial_haplogroup_dict, sample_mapping, edge_mapping, namespace.outfile)
    LOG.info("Finished drawing haplogroups")


def get_arguments() -> argparse.Namespace:
    """Get the arguments provided by the user to this script"""
    parser = argparse.ArgumentParser(description="Erasmus MC: Genetic Identification\n Haplogroup tree creation")

    parser.add_argument("-i", "--input", required=True,
                        help="prediction file (.hg file) containing the predicted haplogroups for different samples. "
                             "Alternalively you can provide a file containing samples in the first column and "
                             "haplogroups in the second.",
                        metavar="FILE")
    parser.add_argument("-c", "--collapse_mode", help="Add this flag to compress the haplogroup tree image and remove"
                                                      " all uninformative haplogroups from it.",
                        action="store_true")
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
    return haplogroups, dict(sample_mapping)


def add_main_haplogroups(haplogroups: List[str]):
    haplogroups.extend([chr(nr) for nr in range(66, 85)])
    haplogroups.append("A0-T")
    haplogroups.append("A00")


def haplogroup_tree_dict(
    haplogroups: List[str]
) -> Dict[str, Set[str]]:
    tree = Tree("Hg_Prediction_tables/tree.json")
    partial_tree_dict = defaultdict(set)
    for name in haplogroups:
        path = []
        node = tree.get(name)
        while node is not None:
            path.append(node.name)
            node = node.parent
        parent = path.pop()
        while len(path) > 0:
            child = path.pop()
            partial_tree_dict[parent].add(child)
            parent = child
    return dict(partial_tree_dict)


def collapse_tree_dict(
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]]
) -> Dict[Tuple[str, str], str]:
    edge_mapping = {}
    remove_keys = set()
    for parent, children in partial_tree_dict.items():
        if parent in remove_keys:
            continue
        orig_parent = parent
        for child in children.copy():
            total_collapsed = 0
            partial_tree_dict[orig_parent].remove(child)
            while True:
                if not can_collapse(child, partial_tree_dict, sample_mapping):
                    break
                remove_keys.add(child)
                child = next(iter(partial_tree_dict[child]))
                total_collapsed += 1
            partial_tree_dict[orig_parent].add(child)
            if total_collapsed > 0:
                edge_mapping[(orig_parent, child)] = str(total_collapsed)
    # remove all the collapsed nodes
    for name in remove_keys:
        del partial_tree_dict[name]
    return edge_mapping


def can_collapse(
    name: str,
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]]
) -> bool:
    # is only a child
    if name not in partial_tree_dict:
        return False
    # more then one child
    if len(partial_tree_dict[name]) != 1:
        return False
    # major node
    if name == 'A0-T' or '-' not in name:
        return False
    # node linked to sample
    if name in sample_mapping:
        return False
    return True


def make_dendrogram(
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]],
    edge_mapping: Dict[Tuple[str, str], str],
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
            edge_name = (parent, child)
            if edge_name not in covered_edges:
                if edge_name in edge_mapping:
                    edge_label = edge_mapping[edge_name]
                else:
                    edge_label = ''
                dot.edge(parent, child, weight="1", label=edge_label)
                covered_edges.add(edge_name)
                if child in sample_mapping:
                    dot.node('\n'.join(sample_mapping[child]), '\n'.join(sample_mapping[child]), shape="box",
                             style='filled', fillcolor="green")
                    dot.edge(child, '\n'.join(sample_mapping[child]), weight="10")
    dot.render(output_file, view=False)
    os.remove(output_file)


if __name__ == '__main__':
    main()
