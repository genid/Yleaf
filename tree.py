# for reading the yfull tree structure
import json


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
