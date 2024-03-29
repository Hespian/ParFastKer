import sys
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib2tikz import save as tikz_save


def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def safe_div(x,y):
    if y == 0:
        return 0
    return x / y

def readGraphFile(inputFile):
    file = open(inputFile, "r")
    graphName = ""
    numEdgesString = ""
    numVerticesString = ""
    iterationStartTimes = []
    iterationStartSizes = []
    partitionFinishTimes = [[] for i in range(32)]
    partitionFinishSizes = [[] for i in range(32)]
    terminationTimes = []
    full_parallel = False
    finish_time = 100000
    for line in file:
        words = line.split()
        if "Number of blocks:" in line:
            full_parallel = int(words[-1]) == 32
        if full_parallel:
            if "Filename:" in line:
                graphName = os.path.os.path.splitext(words[-1])[0]
                print(graphName)
            if "|-Nodes:" in line:
                numVerticesString = words[-1]
            if "|-Edges:" in line:
                numEdgesString = words[-1]
            if "Iteration" in line:
                iterationStartTimes.append(float(words[4]))
                iterationStartSizes.append(int(words[-1]))
            if "Partition" in line:
                partitionNum = int(words[1])
                partitionFinishTimes[partitionNum].append(float(words[6]))
                partitionFinishSizes[partitionNum].append(int(words[-1]))
            if "Termination time" in line:
                terminationTimes.append(float(words[-1]))
            if "Total time spent applying reductions" in line:
                finish_time = float(words[-1])
            if "Kernel size after parallel run" in line:
                iterationStartTimes.append(finish_time)
                iterationStartSizes.append(int(words[-1]))
            if "New repitition: 1" in line:
                break
    file.close()
    plt.yscale('log')
    plt.scatter(iterationStartTimes, iterationStartSizes, color="black", label="Start of iteration", marker="x")
    for time in iterationStartTimes:
        plt.axvline(x=time, color='black', linestyle='--', linewidth=0.5)
    for time in terminationTimes:
        plt.axvline(x=time, color='red', linestyle='--', linewidth=0.5)

    cmap = get_cmap(32)
    cmap = plt.get_cmap('hsv')
    colors = cmap(np.linspace(0, 1.0, 32))
    for i in range(32):
        plt.scatter(partitionFinishTimes[i], partitionFinishSizes[i], color=colors[i], label="partition_" + str(i), marker="x")
    plt.xlim(xmin=0)
    # plt.ylim(ymin=10**4)
    # plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Graph size')
    plt.title(os.path.join(os.path.basename(os.path.split(inputFile)[0]), os.path.basename(inputFile)))
    outputfile = inputFile + '.pdf'
    # plt.show()
    # plt.savefig(outputfile, format='pdf')
    # plt.clf()

def getGraphName(filename):
    res = filename.replace("-sorted", "")
    res = res.replace("-LinearTimeKernel", "")
    res = res.replace(".graph", "")
    return res

directory = sys.argv[1]
outputfile = os.path.join(directory, "plots.pdf")

plt.figure(figsize=(15,30))

graphs = {}
numDirs = 0
for subdir in os.listdir(directory):
    resultsdir = os.path.join(directory, subdir)
    if not os.path.isdir(resultsdir) or resultsdir.endswith("hidden"):
        continue
    numDirs += 1
    for resultsFile in os.listdir(resultsdir):
        resultsFilePath = os.path.join(resultsdir, resultsFile)
        if resultsFile.endswith(".graph") and os.path.isfile(resultsFilePath):
            graphname = getGraphName(resultsFile)
            if not graphname in graphs:
                nextPos = len(graphs) + 1
                graphs[graphname] = nextPos

print(numDirs)
print(len(graphs))
currentdir = 0
for subdir in os.listdir(directory):
    resultsdir = os.path.join(directory, subdir)
    if not os.path.isdir(resultsdir) or resultsdir.endswith("hidden"):
        continue
    for resultsFile in os.listdir(resultsdir):
        resultsFilePath = os.path.join(resultsdir, resultsFile)
        if resultsFile.endswith(".graph") and os.path.isfile(resultsFilePath):
            graphname = getGraphName(resultsFile)
            plt.subplot(len(graphs), numDirs, currentdir + numDirs * (graphs[graphname] - 1) + 1)
            readGraphFile(resultsFilePath)
    currentdir += 1

# plt.subplots_adjust(top=5)
# plt.show()
plt.tight_layout()
plt.savefig(outputfile, format='pdf', height = 100)
