import sys
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib2tikz import save as tikz_save
import get_data_ours



graphs = ["it-2004", "sk-2005", "uk-2007-05", "webbase-2001", "rgg_n26_s0", "del26"]

linearTimeDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/LinearTimeKernels/logs"
partitioningDir = "/home/dhespe/Documents/parallel_reductions/LinearTimeKernels/partitions"
ourTimeDir = "/home/dhespe/Documents/parallel_reductions/results/LinearTimeKernelsScalingAll"

def makePlot(sizes, data, colors, name, big):
    cmapIndex = 1
    for graph in graphs:
        plt.loglog(sizes, data[graph], color=colors[cmapIndex], label=graph, basex=2, marker="x")
        cmapIndex += 1
    plt.axis([1,32,0,10**2])
    plt.xlabel('Number of threads')
    plt.ylabel('Speedup')
    plt.tick_params(axis='y', which='both', left='on', right='on')
    if big:
        plt.legend(loc="upper left")
    plt.rc('font', size=20)
    if sizes[0] != 1:
        sizes.insert(0,1)
    plt.xticks(sizes, sizes)
    plt.savefig(name, bbox_inches="tight")
    plt.clf()

cmap = plt.get_cmap('hsv')
colors = cmap(np.linspace(0, 1.0, len(graphs) + 1))

sizes = []
overall = dict()
quasikernel = dict()
partitioning = dict()
LP = dict()
Rest = dict()
for graph in graphs:
    res = get_data_ours.getOurTimeAndSize(graph, linearTimeDir, partitioningDir, ourTimeDir)
    sizes = sorted(res["scaling_quasikernel_time"].keys())
    overall[graph] = [res["scaling_total"][1] / res["scaling_total"][i] for i in sizes]
    quasikernel[graph] = [res["scaling_quasikernel_time"][1] / res["scaling_quasikernel_time"][i] for i in sizes]
    partitioning[graph] = [res["scaling_partitioning_time"][2] / res["scaling_partitioning_time"][i] for i in sizes[1:]]
    LP[graph] = [res["scaling_LP_time"][1] / res["scaling_LP_time"][i] for i in sizes]
    Rest[graph] = [res["scaling_Rest_time"][1] / res["scaling_Rest_time"][i] for i in sizes]

makePlot(sizes, overall, colors, "Overall.pdf", True)
makePlot(sizes, quasikernel, colors, "Quasikernel.pdf", False)
makePlot(sizes[1:], partitioning, colors, "Partitioning.pdf", False)
makePlot(sizes, LP, colors, "LP.pdf", False)
makePlot(sizes, Rest, colors, "Rest.pdf", False)


