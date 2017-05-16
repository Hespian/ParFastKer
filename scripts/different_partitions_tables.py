import os
import sys
import numpy as np

inputDir = sys.argv[1]
lpaResultsDir = sys.argv[2]

latexFileName = os.path.join(inputDir, "tables.tex")

latexFile = open(latexFileName, 'w')

def getPartitionTime(graph, config):
    filepath = os.path.join(os.path.dirname(graph), "partitions", os.path.basename(graph), config, "16.partition-log")
    file = open(filepath)
    for line in file:
        if ("avg partitioning time elapsed" in line) or ("time spent for partitioning" in line):
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

    for line in file:
        words = line.split()
        if "weight_" in line:
            total_time = 0.0
            partitioned_time = 0.0
            kernel_size = 0
            partition_type = words[0]
        if "Total time spent undoing" in line:
            mean = np.mean(time_per_block)
            time_per_block = []
            if current_rep == num_reps - 1:
                total_time /= num_reps
                partitioned_time /= num_reps
                kernel_size /= num_reps
                if partition_type.count('_') == 2:
                    partitionTime = getPartitionTime(graphResultFile, partition_type)
                    results[partition_type] = (total_time, partitioned_time, round(kernel_size), partitionTime)
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
            time_per_block.append(float(words[-1]))

    file.close()
    return results

def addTotalSummaryTable(results):
    total_time = {}
    partitioned_time = {}
    kernel_size = {}
    for graph_name, graphDict in results.items():
        for config, times in graphDict.items():
            total_time[config] = total_time.get(config, 0.0) + times[0]
            partitioned_time[config] = partitioned_time.get(config, 0.0) + times[1]
            kernel_size[config] = kernel_size.get(config, 0) + times[2]

    base_case = "weight_one_eco"


    latexFile.write("\\begin{tabular}{l|rr}\n")
    latexFile.write("partition type & partitioned time & kernel size \\\\ \n")
    latexFile.write("\\hline \n")

    for config in iter(sorted(total_time.keys())):
        total_time_rel = round(total_time[config] / total_time[base_case] * 100)
        partitioned_time_rel = round(partitioned_time[config] / partitioned_time[base_case] * 100)
        kernel_size_rel = round(kernel_size[config] / kernel_size[base_case] * 100)

        latexFile.write(config.replace("_","\\_") + " & " + str(partitioned_time_rel) + " & " + str(kernel_size_rel) + " \\\\ \n")


    latexFile.write("\\end{tabular} \n")
    latexFile.write("\\vspace{1cm} \n\\newline \n")


def addRelativeSummaryTable(results):
    base_case = "weight_one_eco"
    partitioned_time = {}
    kernel_size = {}
    partition_time = {}
    total_time = {}
    num_graphs = {}
    for graph_name, graphDict in results.items():
        for config, times in graphDict.items():

            partitioned_time_speedup =  (graphDict[base_case][1] / times[1])
            kernel_size_rel = times[2] / graphDict[base_case][2] * 100
            partition_time_speedup = (graphDict[base_case][3] / times[3])
            total_time_speedup = ((graphDict[base_case][1] + graphDict[base_case][3]) / (times[1] + times[3]))

            partitioned_time[config] = partitioned_time.get(config, 0.0) + partitioned_time_speedup
            kernel_size[config] = kernel_size.get(config, 0.0) + kernel_size_rel
            partition_time[config] = partition_time.get(config, 0.0) + partition_time_speedup
            total_time[config] = total_time.get(config, 0.0) + total_time_speedup
            num_graphs[config] = num_graphs.get(config, 0) + 1



    latexFile.write("\\begin{tabular}{ll|rrr|r}\n")
    latexFile.write(" &  & \multicolumn{3}{|c|}{speedup for} & \\\\ \n")
    latexFile.write("weight & configuration & partitioning & kernelization & total & kernel size \\\\ \n")
    latexFile.write("\\hline \n")

    for config in iter(sorted(total_time.keys())):
        partitioned_time_speedup = round(partitioned_time[config] / num_graphs[config], 1)
        kernel_size_rel = round((kernel_size[config] / num_graphs[config]) - 100, 1)
        partition_time_speedup = round(partition_time[config] / num_graphs[config], 1)
        total_time_speedup = round(total_time[config] / num_graphs[config], 1)

        weight = config.split('_')[1].replace("one", "one")
        configuration = config.split('_')[2]

        latexFile.write(weight + " & " + configuration + " & " + str(partition_time_speedup) + " & " + str(partitioned_time_speedup) + " & " + str(total_time_speedup) + " & " + '{0:+}'.format(kernel_size_rel) + "\% \\\\ \n")


    latexFile.write("\\end{tabular} \n")
    latexFile.write("\\vspace{1cm} \n\\newline \n")


def addTable(graphDict):
    base_case = graphDict["weight_one_eco"]

    latexFile.write("\\begin{tabular}{l|rrr}\n")
    latexFile.write("partition type & partitioned time & kernel size\\\\ \n")
    latexFile.write("\\hline \n")

    for partition_type, results in iter(sorted(graphDict.items())):
        total_time_rel = round(results[0] / base_case[0] * 100)
        partitioned_time_rel = round(results[1] / base_case[1] * 100)
        kernel_size_rel = round(results[2] / base_case[2] * 100)

        latexFile.write(partition_type.replace("_","\\_") + " & " + str(partitioned_time_rel) +  " & " + str(kernel_size_rel) + " \\\\ \n")

    latexFile.write("\\end{tabular} \n")
    latexFile.write("\\vspace{1cm} \n\\newline \n")

results = dict()
lpaResults = os.listdir(lpaResultsDir)
for filename in os.listdir(inputDir):
    if filename.endswith(".graph"): 
        results[filename] = readGraphResults(os.path.join(inputDir, filename))
        if filename in lpaResults:
            tempResults = readGraphResults(os.path.join(lpaResultsDir, filename))
            results[filename].update(tempResults)
    else:
        continue

latexFile.write("\\documentclass[varwidth=\\maxdimen, border=10pt]{standalone} \n")
latexFile.write("\\begin{document} \n")

for graph_name, graphDict in results.items():
    latexFile.write(graph_name + "\n")
    latexFile.write("\\newline \n")
    addTable(graphDict)

# latexFile.write("Summary total \n")
# latexFile.write("\\newline \n")
# addTotalSummaryTable(results)

latexFile.write("Summary relative \n")
latexFile.write("\\newline \n")
addRelativeSummaryTable(results)

latexFile.write("\\end{document} \n")
latexFile.close()

os.system("pdflatex -output-directory " + os.path.dirname(latexFileName) + " " + latexFileName)


