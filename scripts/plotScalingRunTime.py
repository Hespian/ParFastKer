import sys
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import *
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib2tikz import save as tikz_save
# from IPython import embed
import get_data_ours
import renameGraphs




graphs = ["it-2004", "uk-2007-05", "webbase-2001",  "del26", "sk-2005", "rgg_n26_s0"]

linearTimeDir = "../../../triangle_counting_paper/MIS_sigmod_pub/results/LinearTimeKernels/logs"
partitioningDir = "../../LinearTimeKernels/partitions"
ourTimeDir = "../../results/LinearTimeKernelsScalingAll"

def makePlot(sizes, data, facecolors, edgecolors, markers, name, ax):
    # if not big:
    cmapIndex = 0
    for graph in graphs:
        ax.loglog(sizes, data[graph], color=edgecolors[cmapIndex], markerfacecolor=facecolors[cmapIndex], markeredgecolor=edgecolors[cmapIndex], label=renameGraphs.renameGraph(graph), basex=2, marker=markers[cmapIndex])
        cmapIndex += 1
    ax.set_xlim([1,32])
    ax.set_ylim([0,10**2])
    ax.set_xlabel('Number of threads', fontsize=14)
    ax.tick_params(axis='y', which='both', left='on', right='on')
    if sizes[0] != 1:
        sizes.insert(0,1)
    ax.set_title(name, fontsize=14)
    ax.set_xticks(sizes)
    ax.set_xticklabels(sizes)

cmap = plt.get_cmap('hsv')
# colors = cmap(np.linspace(0, 1.0, len(graphs) + 1))
facecolors = ["blue", "none", "none", "darkorange", "none", "none"]
edgecolors = ["blue", "green", "red", "darkorange", "magenta", "black"]
markers = ["x", "o", "^", "+" , "s", "d"]

sizes = []
overall = dict()
quasikernel = dict()
partitioning = dict()
LP = dict()
Rest = dict()
for graph in graphs:
    print(graph)
    res = get_data_ours.getOurTimeAndSizeUltrafast(graph, linearTimeDir, partitioningDir, ourTimeDir)
    sizes = sorted(res["scaling_quasikernel_time"].keys())
    overall[graph] = [res["scaling_total"][2] / res["scaling_total"][i] for i in sizes]
    quasikernel[graph] = [res["scaling_quasikernel_time"][2] / res["scaling_quasikernel_time"][i] for i in sizes]
    partitioning[graph] = [res["scaling_partitioning_time"][2] / res["scaling_partitioning_time"][i] for i in sizes[1:]]
    LP[graph] = [res["scaling_LP_time"][2] / res["scaling_LP_time"][i] for i in sizes]
    Rest[graph] = [res["scaling_Rest_time"][2] / res["scaling_Rest_time"][i] for i in sizes]

# print("Overall (32)")
# for graph in graphs:
#     print(overall[graph][-1] / overall[graph][0])

print("----------------------------------")
print("rest (32)")
for graph in graphs:
    print(Rest[graph][-1] / Rest[graph][0])

# print("----------------------------------")
# print("LP (1)")
# for graph in graphs:
#     print(LP[graph][0])


plt.rc('font', size=14)
f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.6, 2.4), sharey=True)
ax1.set_ylabel('Speedup')
makePlot(sizes, overall, facecolors, edgecolors, markers, "Overall", ax1)
# makePlot(sizes, quasikernel, colors, markers, "Quasikernel.pdf", False)
# makePlot(sizes[1:], partitioning, colors, markers, "Partitioning.pdf", False)
makePlot(sizes, Rest, facecolors, edgecolors, markers, "Local Reductions", ax2)
makePlot(sizes, LP, facecolors, edgecolors, markers, "LP Reduction", ax3)
plt.legend(bbox_to_anchor=(-0.7,-0.7), ncol=3, loc='lower center', fontsize=14, frameon=False)
plt.savefig("Scaling.pdf", bbox_inches="tight")
