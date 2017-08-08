import get_data_ours
import get_data_akiba
import get_data_NearLinear
import os
import matplotlib.pyplot as plt

# graphs = ["uk-2002", "arabic-2005", "gsh-2015-tpd", "uk-2005", "it-2004", "sk-2005", "uk-2007-05", "webbase-2001", "asia.osm", "road_usa", "europe.osm", "rgg_n26_s0",  "RHG-100000000-nodes-2000000000-edges", "delaunay_n24", "del26"]
graphs = ["uk-2002", "arabic-2005", "gsh-2015-tpd", "uk-2005", "it-2004", "sk-2005", "uk-2007-05", "webbase-2001", "asia.osm", "road_usa", "europe.osm", "rgg_n26_s0", "delaunay_n24", "del26"]
linearTimeDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/LinearTimeKernels/logs"
nearLinearDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/NearLinear"
partitioningDir = "/home/dhespe/Documents/parallel_reductions/LinearTimeKernels/partitions"
ourTimeDir = "/home/dhespe/Documents/parallel_reductions/results/LinearTimeKernels"
akibaDir = "/home/dhespe/Documents/parallel_reductions/akiba_vertex_cover/results"

def getOurTimeAndSizeSequential(graph):
    res = get_data_ours.getOurTimeAndSize(graph, linearTimeDir, partitioningDir, ourTimeDir)
    result = dict()
    result["time"] = res["sequential_quasikernel_time"] + res["lineartime_time"]
    result["size"] = res["sequential_quasikernel_size"]
    return result

def getOurTimeAndSizeParallel(graph):
    res = get_data_ours.getOurTimeAndSize(graph, linearTimeDir, partitioningDir, ourTimeDir)
    result = dict()
    result["time"] = res["parallel_quasikernel_time"] + res["lineartime_time"] + res["partitioning_time"]
    result["size"] = res["parallel_quasikernel_size"]
    return result

def getAkibaTimeAndSize(graph):
    return get_data_akiba.getAkibaTimeAndSize(graph, akibaDir)

def getNearLinearTimeAndSize(graph):
    return get_data_NearLinear.getNearLinearTimeAndSize(graph, nearLinearDir)

def minProperty(graph, prop):
    oursequential = getOurTimeAndSizeSequential(graph)[prop]
    ourparallel = getOurTimeAndSizeParallel(graph)[prop]
    akiba = getAkibaTimeAndSize(graph)[prop]
    nearLinear = getNearLinearTimeAndSize(graph)[prop]
    data = [oursequential, ourparallel, akiba, nearLinear]
    data = filter(lambda x : x >= 0, data)
    minimum = min(data)
    if minimum == 0:
        return 1
    return minimum

oursizeSequential = []
ourtimeSequential = []
oursizeParallel = []
ourtimeParallel = []
akibasize = []
akibatime = []
nearlinearsize = []
nearlineartime = []

for graph in graphs:
    minsize = minProperty(graph, "size")
    mintime = minProperty(graph, "time")

    oss = getOurTimeAndSizeSequential(graph)["size"] / minsize
    ots = getOurTimeAndSizeSequential(graph)["time"] / mintime
    if oss > 0 and ots > 0:
        oursizeSequential.append(oss)
        ourtimeSequential.append(ots)

    osp = getOurTimeAndSizeParallel(graph)["size"] / minsize
    otp = getOurTimeAndSizeParallel(graph)["time"] / mintime
    if osp > 0 and otp > 0:
        oursizeParallel.append(osp)
        ourtimeParallel.append(otp)

    aks = getAkibaTimeAndSize(graph)["size"] / minsize
    akt = getAkibaTimeAndSize(graph)["time"] / mintime
    if aks > 0 and akt > 0:
        akibasize.append(aks)
        akibatime.append(akt)

    nls = getNearLinearTimeAndSize(graph)["size"] / minsize
    nlt = getNearLinearTimeAndSize(graph)["time"] / mintime
    if nls > 0 and nlt > 0:
        nearlinearsize.append(nls)
        nearlineartime.append(nlt)


# print("We")
# print(oursizeSequential)
# print(ourtimeSequential)

# print("We (parallel)")
# print(oursizeParallel)
# print(ourtimeParallel)

# print("Akiba")
# print(akibasize)
# print(akibatime)

# print("NearLinear")
# print(nearlinearsize)
# print(nearlineartime)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale("log")
ax.set_xscale("log")
ax.scatter(ourtimeSequential, oursizeSequential, label="We (sequential)", marker="x", color="green")
ax.scatter(ourtimeParallel, oursizeParallel, label="We (parallel)", marker="+", color="black")
ax.scatter(akibatime, akibasize, label="Akiba and Iwata", marker="^", edgecolors="blue", facecolors="none")
ax.scatter(nearlineartime, nearlinearsize, label="NearLinear", marker="o", edgecolors="red", facecolors="none")
plt.xlabel("time / fastest time")
plt.ylabel("size / smallest size")
ax.legend()
plt.savefig("summaryplot.pdf")
