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


def main():
    file = "test/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522_chrY/HG00096.mapped.ILLUMINA.bwa.GBR" \
           ".low_coverage.20120522_chrY.out"
    haplotype_dict = read_derived_haplotype_dict(file)
    tree = Tree("Position_files/tree.json")
    get_most_likely_haplotype(tree, haplotype_dict)


def read_derived_haplotype_dict(file):
    # read all the derived lines and count the occurance of each branch
    haplotype_dict = {}
    with open(file) as f:
        f.readline()
        for line in f:
            _, _, _, branches, _, _, _, _, _, _, derived = line.strip().split("\t")
            if derived == "D":

                for branch in branches.split(";"):
                    if branch in haplotype_dict:
                        haplotype_dict[branch] += 1
                    else:
                        haplotype_dict[branch] = 1
    return haplotype_dict


def get_most_likely_haplotype(tree, haplotype_dict):

    sorted_depth_haplotypes = sorted(haplotype_dict.keys(), key=lambda k: tree.get(k).depth, reverse=True)
    not_final_nodes = set()
    scores = {}
    worms = {}
    for haplotype_name in sorted_depth_haplotypes:
        if haplotype_name in not_final_nodes:
            continue
        node = tree.get(haplotype_name)
        parent = node.parent
        score = 0
        worm = f"->{node.name}"
        while parent is not None:
            worm = f"->{parent.name}" + worm
            if parent.name in haplotype_dict:
                score += haplotype_dict[parent.name]
                not_final_nodes.add(parent.name)
            parent = parent.parent
        scores[node.name] = score
        worms[node.name] = worm.replace('->ROOT (Y-Chromosome "Adam")->', "ROOT->")

    best_score = max(scores.items(), key=lambda kv: kv[1])
    print(best_score)
    print(worms[best_score[0]])

    for name in worms:
        print(scores[name])
        print(worms[name])


if __name__ == '__main__':
    main()
