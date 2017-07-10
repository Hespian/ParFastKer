import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save
import math

inputDir = sys.argv[1]

latexFileName = os.path.join(inputDir, "tables.tex")

latexFile = open(latexFileName, 'w')

def getPartitionTime(graph, config):
    filepath = os.path.join(inputDir, "partitions", os.path.basename(graph), config, "16.partition-log")
    file = open(filepath)
    for line in file:
        if "avg partitioning time elapsed" in line:
            words = line.split()
            return float(words[-1])

# total_time, partitioned_time, kernel_size
def readGraphResults(graphResultFile):
    file = open(graphResultFile, "r")
    partition_type = ""
    num_reps = 0;
    total_time = 0.0
    partitioned_time = 0.0;
    graphName = ""
    current_rep = 0
    kernel_size = 0
    stddev = 0.0

    time_per_block = []

    results = dict()
    max_num_blocks = 0

    for line in file:
        words = line.split()
        if "Number of blocks:" in line:
            num_blocks = int(words[-1])
            time_per_block = [0.0] * num_blocks
        if "weight_" in line:
            total_time = 0.0
            partitioned_time = 0.0
            kernel_size = 0
            partition_type = words[0]
        if "Total time spent undoing" in line:
            if current_rep == num_reps - 1:
                time_per_block = [x / num_reps for x in time_per_block]
                total_time /= num_reps
                partitioned_time /= num_reps
                kernel_size /= num_reps
                if partition_type.count('_') == 2:
                    # partitionTime = getPartitionTime(graphResultFile, partition_type)
                    if max_num_blocks < num_blocks:
                        max_num_blocks = num_blocks
                        results[partition_type] = (time_per_block)
        if "Filename:" in line:
            graphName = words[-1]
        if "Number of repititions" in line:
            num_reps = int(words[-1])
        if "New repitition:" in line:
            current_rep = int(words[-1])
        if "Total time spent applying reductions  : " in line:
            total_time += float(words[6])
        if "Kernel size after parallel run: " in line:
            kernel_size += int(words[5])
        if "Rest time" in line:
            partitioned_time += float(words[-1])
        if "Time spent applying reductions" in line:
            block = int(words[0][:-1])
            time_per_block[block] += (float(words[-1]))

    file.close()
    return results

fig = plt.figure(figsize=(10,8))
current_graph = 1
def drawBoxplots(graphResults, graphname, inputdir):
    graphs = list(results.keys())
    data =  [graphResults["weight_one_ultrafast"]]
    labels =  ["uniform"]

    global current_graph
    dimensions = math.ceil(math.sqrt(numGraphs))
    ax = fig.add_subplot(dimensions,dimensions,current_graph)
    current_graph = current_graph + 1
    plt.boxplot(data, widths = 0.7, showmeans=True)
    # plt.xlabel('Weight')
    # plt.ylabel('Run time (s)')
    plt.title(graphname)
    # ax.set_xticklabels(labels=labels, rotation=45)
    # plt.savefig(os.path.join(plotsDir, "runtime_boxplots"))
    # plt.clf()

    filename = os.path.join(inputdir, 'boxplot' + graphname)
    # tikz_save(filename + '.tikz', figureheight = '\\figureheight', figurewidth = '\\figurewidth')
    # plt.show()
    # plt.savefig(filename+ ".boxplot.pdf")
    # plt.clf()

results = dict()
for filename in os.listdir(inputDir):
    if filename.endswith(".graph"): 
        results[filename] = readGraphResults(os.path.join(inputDir, filename))
        continue
    else:
        continue

numGraphs = len(results)
for name, graph in sorted(results.items()):
    drawBoxplots(graph, name, inputDir)

# plt.show()
plt.savefig(os.path.join(inputDir, "boxplots"))


