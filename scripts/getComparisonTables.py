import get_data_ours
import get_data_akiba
import get_data_NearLinear
from tabulate import tabulate
import os
import csv


graphs = ["uk-2002", "arabic-2005", "gsh-2015-tpd", "uk-2005", "it-2004", "sk-2005", "uk-2007-05", "webbase-2001", "asia.osm", "road_usa", "europe.osm", "rgg_n26_s0",  "RHG-100000000-nodes-2000000000-edges", "delaunay_n24", "del26"]
linearTimeDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/LinearTimeKernels/logs"
nearLinearDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/NearLinear"
partitioningDir = "/home/dhespe/Documents/parallel_reductions/LinearTimeKernels/partitions"
ourTimeDir = "/home/dhespe/Documents/parallel_reductions/results/LinearTimeKernels"
akibaDir = "/home/dhespe/Documents/parallel_reductions/akiba_vertex_cover/results"

def getOurTimeForSize(graph, targetSize):
    return get_data_ours.getOurTimeForSize(graph, linearTimeDir, partitioningDir, ourTimeDir, targetSize)

def getOurTimeAndSize(graph):
    return get_data_ours.getOurTimeAndSize(graph, linearTimeDir, partitioningDir, ourTimeDir)

def getAkibaTimeAndSize(graph):
    return get_data_akiba.getAkibaTimeAndSize(graph, akibaDir)

def getAkibaTimeForSize(graph, targetsize):
    return get_data_akiba.getAkibaTimeForSize(graph, akibaDir, targetsize)

def getNearLinearTimeAndSize(graph):
    return get_data_NearLinear.getNearLinearTimeAndSize(graph, nearLinearDir)

def RemoveNegatives(line):
    result = []
    for item in line:
        stringitem = ""
        if isinstance(item, int):
            if item < 0:
                stringitem = "*"
            else:
                stringitem = "\\numprint{" + str(item) + "}"
        elif isinstance(item, float):
            if item <= 0:
                stringitem = "*"
            else:
                stringitem = "\\numprint{%.3f}" % item
        elif isinstance(item, str):
            stringitem = item.replace("_", "\\_")

        if item ==   "RHG-100000000-nodes-2000000000-edges":
            stringitem = "RHG-100M-2G"
        result.append(stringitem)
    return result

def writeToFile(headers, data, filename):
    ofile = open(filename, "w")
    writer = csv.writer(ofile)
    writer.writerow(headers)
    for line in data:
        writer.writerow(line)
    ofile.close()

def ourAlgorithmSequential():
    data = []
    headers = ["graph", "ourSize", "LinearTimeTime", "quasiKernelTime", "quasiKernelTotal", "kernelSize", "kernelTime", "kerneltotal"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["sequential_quasikernel_time"]
        total_time_kernel = total_time + ourTimeAndSize["sequential_kernel_time"]
        line = [graph, int(ourTimeAndSize["sequential_quasikernel_size"]), ourTimeAndSize["lineartime_time"], ourTimeAndSize["sequential_quasikernel_time"], total_time, int(ourTimeAndSize["sequential_kernel_size"]), ourTimeAndSize["sequential_kernel_time"], total_time_kernel]
        data.append(RemoveNegatives(line))

    writeToFile(headers, data, "ourAlgorithmSequential.csv")

def ourAlgorithmParallel():
    data = []
    headers = ["graph", "ourSize", "LinearTimeTime", "partitioningTime", "quasiKernelTime", "quasiKernelTotal", "kernelSize", "kernelTime", "kerneltotal"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["partitioning_time"] + ourTimeAndSize["parallel_quasikernel_time"]
        total_time_kernel = total_time + ourTimeAndSize["parallel_kernel_time"]
        line = [graph, int(ourTimeAndSize["parallel_quasikernel_size"]), ourTimeAndSize["lineartime_time"], ourTimeAndSize["partitioning_time"], ourTimeAndSize["parallel_quasikernel_time"], total_time, int(ourTimeAndSize["parallel_kernel_size"]), ourTimeAndSize["parallel_kernel_time"], total_time_kernel]
        data.append(RemoveNegatives(line))

    writeToFile(headers, data, "ourAlgorithmParallel.csv")

def sequentialQuasikernelComparsionAkiba():
    data = []
    headers = ["graph", "AkibaSize", "AkibaTime", "ourSize", "LinearTimeTime", "ourTime", "quasiKernelTotal", "speedupTotal", "sameSizeTime", "speedupSameSize"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        akibaTimeAndSize = getAkibaTimeAndSize(graph)
        akibaTimeForOurSize = getAkibaTimeForSize(graph, ourTimeAndSize["sequential_quasikernel_size"])
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["sequential_quasikernel_time"]
        speedup_total = akibaTimeAndSize["time"] / total_time
        same_size_speedup = akibaTimeForOurSize / total_time
        line = [graph, akibaTimeAndSize["size"], akibaTimeAndSize["time"], int(ourTimeAndSize["sequential_quasikernel_size"]), ourTimeAndSize["lineartime_time"], ourTimeAndSize["sequential_quasikernel_time"], total_time, speedup_total, akibaTimeForOurSize, same_size_speedup]
        data.append(RemoveNegatives(line))

    writeToFile(headers, data, "sequentialQuasiKernelComparisonAkiba.csv")
    # return(tabulate(data, headers=headers, floatfmt=".3f", tablefmt="latex"))


def sequentialKernelComparsionAkiba():
    data = []
    headers = ["graph", "Akiba size", "Akiba time", "our size", "LinearTime time", "our time quasikernel", "our time kernel", "total", "speedup", "same size time", "speedup"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        akibaTimeAndSize = getAkibaTimeAndSize(graph)
        akibaTimeForOurSize = getAkibaTimeForSize(graph, ourTimeAndSize["sequential_kernel_size"])
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["sequential_quasikernel_time"] + ourTimeAndSize["sequential_kernel_time"]
        speedup_total = akibaTimeAndSize["time"] / total_time
        same_size_speedup = akibaTimeForOurSize / total_time
        line = [graph, akibaTimeAndSize["size"], akibaTimeAndSize["time"], int(ourTimeAndSize["sequential_kernel_size"]), ourTimeAndSize["lineartime_time"], ourTimeAndSize["sequential_quasikernel_time"], ourTimeAndSize["sequential_kernel_time"],total_time, speedup_total, akibaTimeForOurSize, same_size_speedup]
        data.append(line)

    # return(tabulate(data, headers=headers, floatfmt=".3f", tablefmt="latex"))


def parallelKernelComparsionAkiba():
    data = []
    headers = ["graph", "Akiba size", "Akiba time", "our size", "LinearTime time", "partitioning time", "our time quasikernel", "our time kernel", "total", "speedup", "same size time", "speedup"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        akibaTimeAndSize = getAkibaTimeAndSize(graph)
        akibaTimeForOurSize = getAkibaTimeForSize(graph, ourTimeAndSize["parallel_kernel_size"])
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["partitioning_time"]+ ourTimeAndSize["parallel_quasikernel_time"] + ourTimeAndSize["parallel_kernel_time"]
        speedup_total = akibaTimeAndSize["time"] / total_time
        same_size_speedup = akibaTimeForOurSize / total_time
        line = [graph, akibaTimeAndSize["size"], akibaTimeAndSize["time"], int(ourTimeAndSize["parallel_kernel_size"]), ourTimeAndSize["lineartime_time"], ourTimeAndSize["partitioning_time"],ourTimeAndSize["parallel_quasikernel_time"], ourTimeAndSize["parallel_kernel_time"],total_time, speedup_total, akibaTimeForOurSize, same_size_speedup]
        data.append(line)

    # return(tabulate(data, headers=headers, floatfmt=".3f", tablefmt="latex"))

def parallelQuasikernelComparsionAkiba():
    data = []
    headers = ["graph", "AkibaSize", "AkibaTime", "ourSize", "ourTimeQuasiKernelTotal", "totalSpeedup", "sameSizeTime", "sameSizeSpeedup"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        akibaTimeAndSize = getAkibaTimeAndSize(graph)
        akibaTimeForOurSize = getAkibaTimeForSize(graph, ourTimeAndSize["parallel_quasikernel_size"])
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["partitioning_time"]+ ourTimeAndSize["parallel_quasikernel_time"]
        speedup_total = akibaTimeAndSize["time"] / total_time
        same_size_speedup = akibaTimeForOurSize / total_time
        line = [graph, akibaTimeAndSize["size"], akibaTimeAndSize["time"], int(ourTimeAndSize["parallel_quasikernel_size"]), total_time, speedup_total, akibaTimeForOurSize, same_size_speedup]
        data.append(RemoveNegatives(line))

    writeToFile(headers, data, "parallelQuasiKernelComparisonAkiba.csv")

    # return(tabulate(data, headers=headers, floatfmt=".3f", tablefmt="latex"))


def sequentialQuasikernelComparsionNearLinear():
    data = []
    headers = ["graph", "NearLinearSize", "NearLinearTime", "ourSize", "ourTimetotal", "totalSpeedup", "sameSizeTotalTime", "sameSizeSpeedup"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        nearLinearTimeAndSize = getNearLinearTimeAndSize(graph)
        ourTimeForNearLinearSize = getOurTimeForSize(graph, nearLinearTimeAndSize["size"])["sequential"]
        totalTimeSameSize = ourTimeForNearLinearSize + ourTimeAndSize["lineartime_time"]
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["sequential_quasikernel_time"]
        speedup_total = nearLinearTimeAndSize["time"] / total_time
        same_size_speedup = nearLinearTimeAndSize["time"] / totalTimeSameSize
        if totalTimeSameSize < 0:
            same_size_speedup = -1
        line = [graph, nearLinearTimeAndSize["size"], nearLinearTimeAndSize["time"], int(ourTimeAndSize["sequential_quasikernel_size"]), total_time, speedup_total, totalTimeSameSize, same_size_speedup]
        data.append(RemoveNegatives(line))

    writeToFile(headers, data, "sequentialQuasiKernelComparisonNearLinear.csv")

    # return(tabulate(data, headers=headers, floatfmt=".3f", tablefmt="latex"))


def parallelQuasikernelComparsionNearLinear():
    data = []
    headers = ["graph", "NearLinearSize", "NearLinearTime", "ourSize", "ourTimetotal", "totalSpeedup", "sameSizeTotalTime", "sameSizeSpeedup"]
    for graph in graphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        nearLinearTimeAndSize = getNearLinearTimeAndSize(graph)
        ourTimeForNearLinearSize = getOurTimeForSize(graph, nearLinearTimeAndSize["size"])["parallel"]
        totalTimeSameSize = ourTimeForNearLinearSize + ourTimeAndSize["partitioning_time"] + ourTimeAndSize["lineartime_time"]
        total_time = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["partitioning_time"] + ourTimeAndSize["parallel_quasikernel_time"]
        speedup_total = nearLinearTimeAndSize["time"] / total_time
        same_size_speedup = nearLinearTimeAndSize["time"] / totalTimeSameSize
        if totalTimeSameSize < 0:
            same_size_speedup = -1
        line = [graph, nearLinearTimeAndSize["size"], nearLinearTimeAndSize["time"], int(ourTimeAndSize["parallel_quasikernel_size"]), total_time, speedup_total, totalTimeSameSize, same_size_speedup]
        data.append(RemoveNegatives(line))

    writeToFile(headers, data, "parallelQuasiKernelComparisonNearLinear.csv")
    # return(tabulate(data, headers=headers, floatfmt=".3f", tablefmt="latex"))

def writeTableToLatexFile(table, filename):
    file = open(filename, "w")
    file.write("\\documentclass[varwidth=\\maxdimen, border=10pt]{standalone} \n")
    file.write("\\begin{document} \n")
    file.write(table)
    file.write("\\end{document} \n")
    file.close()
    os.system("pdflatex " + filename)

# writeTableToLatexFile(sequentialQuasikernelComparsionAkiba(), "sequentialQuasikernelComparsionAkiba.tex")
# writeTableToLatexFile(sequentialKernelComparsionAkiba(), "sequentialKernelComparsionAkiba.tex")
# writeTableToLatexFile(parallelKernelComparsionAkiba(), "parallelKernelComparsionAkiba.tex")
# writeTableToLatexFile(parallelQuasikernelComparsionAkiba(), "parallelQuasikernelComparsionAkiba.tex")
# writeTableToLatexFile(sequentialQuasikernelComparsionNeaLinear(), "sequentialQuasikernelComparsionNeaLinear.tex")
# writeTableToLatexFile(parallelQuasikernelComparsionNeaLinear(), "parallelQuasikernelComparsionNeaLinear.tex")
sequentialQuasikernelComparsionAkiba()
ourAlgorithmSequential()
ourAlgorithmParallel()
parallelQuasikernelComparsionAkiba()
parallelQuasikernelComparsionNearLinear()
sequentialQuasikernelComparsionNearLinear()
