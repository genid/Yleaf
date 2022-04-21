#!/usr/bin/env python

"""
Helper functions for reading a Json file containing the YFull tree structure

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""


import json
from typing import TYPE_CHECKING, Dict, List, Union

if TYPE_CHECKING:
    from pathlib import Path


class Tree:
    """Hold the YFull tree in a quickly acccesible structure"""
    __ROOT_KEY: str = 'ROOT (Y-Chromosome "Adam")'

    node_mapping: Dict[str, "Node"]

    def __init__(
        self,
        file: Union["Path", str]
    ):
        self.node_mapping = {}
        self._construct_tree(file)

    def _construct_tree(
        self,
        file: Union["Path", str]
    ):
        """Construct the tree based on a Json file containing a dictionary keyed on node names with lists of child
        nodes as values"""
        with open(file) as f:
            tree_dict = json.load(f)

        # base has no parent
        base_node = Node(self.__ROOT_KEY, None, 0, ["A00", "A0-T"])
        self.node_mapping[self.__ROOT_KEY] = base_node
        self._recursive_read(tree_dict, base_node, 0)

    def _recursive_read(
        self,
        tree_dict: Dict[str, List[str]],
        node: "Node", depth: int
    ):
        """Function for recursively reading a dictionary keyed on node names with child nodes as values"""
        depth += 1
        for node_name in tree_dict[node.name]:
            children = []
            if node_name in tree_dict:
                children = tree_dict[node_name]
            new_node = Node(node_name, node, depth, children)
            self.node_mapping[node_name] = new_node
            self._recursive_read(tree_dict, new_node, depth)

    def get(
        self,
        node_name: str
    ) -> "Node":
        """Get a node from the tree based on string name"""
        return self.node_mapping[node_name]


class Node:
    """Nodes used in the tree in order to make accesing information a bit more clear as well as protect the
    modification of variables trough properties"""
    _name: str
    _parent: Union["Node", None]
    _depth: int
    _children: List[str]

    __slots__ = "_name", "_parent", "_depth", "_children"

    def __init__(
        self,
        name: str,
        parent_node: Union["Node", None],
        depth: int,
        children: List[str]
    ):
        self._name = name
        self._parent = parent_node
        self._depth = depth
        self._children = children

    @property
    def name(self) -> str:
        return self._name

    @property
    def parent(self) -> Union["Node", None]:
        return self._parent

    @property
    def depth(self) -> int:
        return self._depth

    @property
    def children(self) -> List[str]:
        return self._children
