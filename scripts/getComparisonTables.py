import get_data_ours
import get_data_akiba
import get_data_NearLinear
import get_data_LinearTime
import renameGraphs
from tabulate import tabulate
import os
import csv


graphs = ["uk-2002", "arabic-2005", "gsh-2015-tpd", "uk-2005", "it-2004", "sk-2005", "uk-2007-05", "webbase-2001", "asia.osm", "road_usa", "europe.osm", "rgg_n26_s0",  "RHG-100000000-nodes-2000000000-edges", "delaunay_n24", "del26"]
summaryTableGraphs = ["it-2004", "webbase-2001", "asia.osm", "europe.osm", "rgg_n26_s0", "RHG-100000000-nodes-2000000000-edges", "del26"]
linearTimeDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/LinearTimeKernels/logs"
nearLinearDir = "/home/dhespe/Documents/triangle_counting_paper/MIS_sigmod_pub/results/NearLinear"
partitioningDir = "/home/dhespe/Documents/parallel_reductions/LinearTimeKernels/partitions"
ourTimeDir = "/home/dhespe/Documents/parallel_reductions/results/LinearTimeKernels"
akibaDir = "/home/dhespe/Documents/parallel_reductions/akiba_vertex_cover/results"
graphdir = "/home/dhespe/Documents/parallel_reductions/graphs"
noReductionStoppingDir = "/home/dhespe/Documents/parallel_reductions/results/LinearTimeKernelsNoReductionStopping"

def getGraphSize(graph):
    result = dict()
    for filename in os.listdir(graphdir):
        if filename.endswith(".graph") and filename.startswith(graph):
            with open(os.path.join(graphdir, filename), "r") as file:
                for line in file:
                    if not line.startswith("%"):
                        words = line.split()
                        return int(words[0])

def getOurTimeForSize(graph, targetSize):
    return get_data_ours.getOurTimeForSize(graph, linearTimeDir, partitioningDir, ourTimeDir, targetSize)

def getOurTimeAndSize(graph):
    return get_data_ours.getOurTimeAndSizeUltrafast(graph, linearTimeDir, partitioningDir, ourTimeDir)

def getAkibaTimeAndSize(graph):
    return get_data_akiba.getAkibaTimeAndSize(graph, akibaDir)

def getAkibaTimeForSize(graph, targetsize):
    return get_data_akiba.getAkibaTimeForSize(graph, akibaDir, targetsize)

def getNearLinearTimeAndSize(graph):
    return get_data_NearLinear.getNearLinearTimeAndSize(graph, nearLinearDir)

def getLinearTimeTimeAndSize(graph):
    return get_data_LinearTime.getLinearTimeTimeAndSize(graph, linearTimeDir)

def RemoveNegatives(line, decimalPlaces=1):
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
                formatstring = "\\numprint{%." + str(decimalPlaces) + "f}"
                stringitem = formatstring % item
        elif isinstance(item, str):
            stringitem = renameGraphs.renameGraph(item)
            stringitem = stringitem.replace("_", "\\_")
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

def percentOf(value, of):
    return value * 100 / of

def scaleToSameOrder(sizes):
    minSize = min(filter(lambda x: x >= 0, sizes))
    if minSize >= 100000:
        return ["\\numprint{%.1f}M" % (size / 1000000) if size >= 0 else "*" for size in sizes]
    elif minSize >= 1000:
        return ["\\numprint{%.1f}K" % (size / 1000) if size >= 0 else "*" for size in sizes]
    else:
        return ["\\numprint{"+ str(int(size)) + "}" if size >= 0 else "*" for size in sizes]

def addNumprint(numbers):
        return ["\\numprint{%.1f}" % number if number >= 0 else "*" for number in numbers]

def markSmallest(numbers, numstrings, textstyle):
    minNumber = min(filter(lambda x: x >= 0, numbers))
    flags = [x == minNumber for x in numbers]
    return[textstyle + "{" + numstring + "}" if flag else numstring for flag, numstring in zip(flags, numstrings)]


def getSummaryComparison(graphsForTable, filename):
    data = []
    headers = ["graph", "graphSize", "LinearTimeTime", "LinearTimeSize", "NearLinearTime", "NearLinearSize", "AkibaTime", "AkibaSize", "OurSequentialTime", "OurSequentialSize", "OurParallelTime", "OurParallelSize", "ParallelSpeedupAkiba"]

    for graph in graphsForTable:
        ourTimeAndSize = getOurTimeAndSize(graph)
        nearLinearTimeAndSize = getNearLinearTimeAndSize(graph)
        linearTimeTimeAndSize = getLinearTimeTimeAndSize(graph)
        akibaTimeAndSize = getAkibaTimeAndSize(graph)
        total_time_parallel = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["partitioning_time"] + ourTimeAndSize["parallel_quasikernel_time"]
        total_time_sequential = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["sequential_quasikernel_time"]
        parallel_speedup_over_akiba = akibaTimeAndSize["time"] / total_time_parallel
        graphSize = getGraphSize(graph)
        graphSizeString = "\\numprint{" + str(int(round(graphSize / 1000000, 0))) + "}M"

        linearTimeSizeString, nearLinearSizeString, AkibaSizeString, SequentialSizeString, ParallelSizeString = scaleToSameOrder([linearTimeTimeAndSize["size"], nearLinearTimeAndSize["size"], akibaTimeAndSize["size"], ourTimeAndSize["sequential_quasikernel_size"], ourTimeAndSize["parallel_quasikernel_size"]])
        linearTimeSizeString, nearLinearSizeString, AkibaSizeString, SequentialSizeString, ParallelSizeString = markSmallest([linearTimeTimeAndSize["size"], nearLinearTimeAndSize["size"], akibaTimeAndSize["size"], ourTimeAndSize["sequential_quasikernel_size"], ourTimeAndSize["parallel_quasikernel_size"]], [linearTimeSizeString, nearLinearSizeString, AkibaSizeString, SequentialSizeString, ParallelSizeString], "\\textit")

        linearTimeTimeString, nearLinearTimeString, AkibaTimeString, SequentialTimeString, ParallelTimeString = addNumprint([linearTimeTimeAndSize["time"], nearLinearTimeAndSize["time"], akibaTimeAndSize["time"], total_time_sequential, total_time_parallel])
        linearTimeTimeString, nearLinearTimeString, AkibaTimeString, SequentialTimeString, ParallelTimeString = markSmallest([linearTimeTimeAndSize["time"], nearLinearTimeAndSize["time"], akibaTimeAndSize["time"], total_time_sequential, total_time_parallel], [linearTimeTimeString, nearLinearTimeString, AkibaTimeString, SequentialTimeString, ParallelTimeString], "\\textbf")

        line = [graph, graphSizeString, linearTimeTimeString, linearTimeSizeString, nearLinearTimeString, nearLinearSizeString, AkibaTimeString, AkibaSizeString, SequentialTimeString, SequentialSizeString, ParallelTimeString, ParallelSizeString, parallel_speedup_over_akiba]
        data.append(RemoveNegatives(line, 1))

    writeToFile(headers, data, filename)

def getSummaryComparisonPercent():
    data = []
    headers = ["graph", "graphSize", "NearLinearTime", "NearLinearSize", "AkibaTime", "AkibaSize", "OurSequentialTime", "OurSequentialSize", "OurParallelTime", "OurParallelSize", "ParallelSpeedupAkiba"]

    for graph in summaryTableGraphs:
        ourTimeAndSize = getOurTimeAndSize(graph)
        nearLinearTimeAndSize = getNearLinearTimeAndSize(graph)
        akibaTimeAndSize = getAkibaTimeAndSize(graph)
        total_time_parallel = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["partitioning_time"] + ourTimeAndSize["parallel_quasikernel_time"]
        total_time_sequential = ourTimeAndSize["lineartime_time"] + ourTimeAndSize["sequential_quasikernel_time"]
        parallel_speedup_over_akiba = akibaTimeAndSize["time"] / total_time_parallel
        graphSize = getGraphSize(graph)
        graphSizeString = "\\numprint{%.1f} M" % (graphSize / 1000000)
        nearLinearSizeString = percentOf(nearLinearTimeAndSize["size"], graphSize)
        AkibaSizeString = percentOf(akibaTimeAndSize["size"], graphSize)
        SequentialSizeString = percentOf(ourTimeAndSize["sequential_quasikernel_size"], graphSize)
        ParallelSizeString = percentOf(ourTimeAndSize["parallel_quasikernel_size"], graphSize)
        line = [graph, graphSizeString, nearLinearTimeAndSize["time"], nearLinearSizeString, akibaTimeAndSize["time"], AkibaSizeString, total_time_sequential, SequentialSizeString, total_time_parallel, ParallelSizeString, parallel_speedup_over_akiba]
        data.append(RemoveNegatives(line, 2))

    writeToFile(headers, data, "summaryComparisonPercent.csv")


def reductionStoppingcomparison():
    data = []
    headers = ["graph", "sizeDiff", "speedup"]

    for graph in graphs:
        resultsWith = getOurTimeAndSize(graph)
        resultsWithout = get_data_ours.getOurTimeAndSizeUltrafast(graph, linearTimeDir, partitioningDir, noReductionStoppingDir)
        timeWith = resultsWith["lineartime_time"] + resultsWith["partitioning_time"] + resultsWith["parallel_quasikernel_time"]
        timeWithout = resultsWithout["lineartime_time"] + resultsWithout["partitioning_time"] + resultsWithout["parallel_quasikernel_time"]
        speedup = timeWithout / timeWith
        sizeWith = resultsWith["parallel_quasikernel_size"]
        sizeWithout = resultsWithout["parallel_quasikernel_size"]
        sizeDifference = sizeWithout - sizeWith
        sizeDifference /= float(sizeWithout)
        sizeDifference *= 100
        if round(speedup, 1) == 1.0:
            continue
        sizeDiffString = "\\numprint{%.1f}" % sizeDifference
        line = [graph, sizeDiffString, speedup]
        data.append(RemoveNegatives(line,1))
        # data.append(line)

    writeToFile(headers, data, "reductionStoppingComparison.csv")
    # print(tabulate(data, headers=headers, floatfmt=".3f"))

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
getSummaryComparison(summaryTableGraphs ,"summaryComparison.csv")
getSummaryComparison(graphs ,"allGraphsComparison.csv")
getSummaryComparisonPercent()
reductionStoppingcomparison()

# for graph in graphs:
#     ourTimeAndSize = getOurTimeAndSize(graph)
#     akibaTimeAndSize = getAkibaTimeAndSize(graph)
#     graphSize = getGraphSize(graph)
#     sizediff = (ourTimeAndSize["parallel_quasikernel_size"] - akibaTimeAndSize["size"]) / graphSize
#     print(graph + " " + str(sizediff))
