import sys
import os
from tabulate import tabulate

def get_size(graph_file):
    graph = open(graph_file, "r")
    for line in graph:
        if line[0] != "%":
            # print(line.split())
            return int(line.split()[0]), int(line.split()[1])

graphdir = sys.argv[1]
print("graphdir = "+ graphdir)
results = []
for filename in os.listdir(graphdir):
    graph_file = os.path.join(graphdir, filename)
    if not os.path.isfile(graph_file) or not graph_file.endswith(".graph"):
        continue
    nodes, edges = get_size(graph_file)
    if nodes >= 10000000:
        # results.append([ filename, '{0:,}'.format(nodes), '{0:,}'.format(edges) ])
        results.append([ filename, (nodes), (edges) ])

print( tabulate(results, headers=["graph", "nodes", "edges"], tablefmt="latex") )
