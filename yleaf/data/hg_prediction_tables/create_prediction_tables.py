import networkx as nx
import json

major_tree = json.load(open('tree_major.json'))

G = nx.Graph(major_tree)

# Get all the leaves
leaves = [n for n in G.nodes()]

print(leaves)

# for leaf in leaves:
#     # Get the leaves from the path from the leaf to the root
#     path = nx.shortest_path(G, leaf, "ROOT (Y-Chromosome \"Adam\")")
#     der_list = [n for n in path if n != "ROOT (Y-Chromosome \"Adam\")"]
#     anc_list = [n for n in leaves if n not in der_list and n != "ROOT (Y-Chromosome \"Adam\")"]
#
#     with open(f'major_tables/{leaf}_int.txt', 'w') as f:
#         for der in der_list:
#             f.write(f'{der}\tD\n')
#         for anc in anc_list:
#             f.write(f'{anc}\tA\n')
