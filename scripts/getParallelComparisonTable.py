import sys
import os

def getPartitionTime(graph, numBlocks):
    filepath = os.path.join(os.path.dirname(graph), "partitions", os.path.basename(graph), "weight_one_ultrafast", str(numBlocks) + ".partition-log")
    file = open(filepath)
    for line in file:
        if ("avg partitioning time elapsed" in line) or ("time spent for partitioning" in line):
            words = line.split()
            return float(words[-1])
    print("Not time found in " + filepath)

def getAkibaRuntimeForSizeFromFile(filename, targetSize):
    file = open(filename)
    for line in file:
        if "Current size" in line:
            words = line.split()
            size = int(words[-1])
            if size <= targetSize:
                return float(words[0][:-1]) / 1000.0
    print(filename)
    raise

def getAkibaRuntimeForSize(akibaResultsDir, graphname, targetSize):
    numFiles = 0
    time = 0
    for filename in os.listdir(akibaResultsDir):
        if graphname in filename:
            numFiles += 1
            time += getAkibaRuntimeForSizeFromFile(os.path.join(akibaResultsDir, filename), targetSize)
    return time / numFiles

def getAkibaFinalTimeAndSizeFromFile(filename):
    file = open(filename)
    for line in file:
        if "Final kernel size:" in line:
            words = line.split()
            size = int(words[-1])
            time = float(words[0][:-1]) / 1000.0
            return (time, size)
    print(filename)
    raise

def getAkibaFinalTimeAndSize(akibaResultsDir, graphname):
    numFiles = 0
    time = 0
    size = 0
    for filename in os.listdir(akibaResultsDir):
        if graphname in filename:
            numFiles += 1
            temptime, tempsize = getAkibaFinalTimeAndSizeFromFile(os.path.join(akibaResultsDir, filename))
            time += temptime
            size += tempsize
    return (time / numFiles, size / numFiles)


def getOurTimeAndSize(inputFile):
    file = open(inputFile, "r")
    sizes = []
    parallel = True
    block_sizes = {}
    time_per_block = {}
    parallel_time = {}
    parallel_size = {}
    num_isolated_clique_reductions = {}
    num_vertex_fold_reductions = {}
    num_twins_removed = {}
    num_twins_folded = {}
    num_unconfined_reductions = {}
    num_lp_removed = {}
    num_diamond_reductions = {}
    sequential_time = {}
    num_blocks = 0
    num_reps = 0
    current_rep = 0
    graphName = ""
    numEdgesString = ""
    numVerticesString = ""
    for line in file:
        words = line.split()
        if "Total time spent undoing" in line:
            if current_rep == num_reps - 1:
                for block_num in range(num_blocks):
                    time_per_block[num_blocks][block_num] /= num_reps
                    num_isolated_clique_reductions[num_blocks][block_num] /= num_reps
                    num_vertex_fold_reductions[num_blocks][block_num] /= num_reps
                    num_twins_removed[num_blocks][block_num] /= num_reps
                    num_twins_folded[num_blocks][block_num] /= num_reps
                    num_unconfined_reductions[num_blocks][block_num] /= num_reps
                    num_diamond_reductions[num_blocks][block_num] /= num_reps
                num_lp_removed[num_blocks] /= num_reps
                parallel_time[num_blocks] /= num_reps
                parallel_size[num_blocks] /= num_reps
                sequential_time[num_blocks] /= num_reps
        if "Filename:" in line:
            graphName = os.path.os.path.splitext(words[-1])[0]
            print(graphName)
        if "|-Nodes:" in line:
            numVerticesString = words[-1]
        if "|-Edges:" in line:
            numEdgesString = words[-1]
        if "Number of repititions" in line:
            num_reps = int(words[-1])
        if "Number of blocks:" in line:
            num_blocks = int(words[-1])
            sizes.append(num_blocks)
            block_sizes[num_blocks] = [0] * num_blocks
            time_per_block[num_blocks] = [0.0] * num_blocks
            num_isolated_clique_reductions[num_blocks] = [0] * num_blocks
            num_vertex_fold_reductions[num_blocks] = [0] * num_blocks
            num_twins_removed[num_blocks] = [0] * num_blocks
            num_twins_folded[num_blocks] = [0] * num_blocks
            num_unconfined_reductions[num_blocks] = [0] * num_blocks
            num_diamond_reductions[num_blocks] = [0] * num_blocks
            parallel_time[num_blocks] = 0
            parallel_size[num_blocks] = 0
            sequential_time[num_blocks] = 0
            num_lp_removed[num_blocks] = 0
        if "New repitition:" in line:
            parallel = True
            current_rep = int(words[-1])
        if "Before call to sequential reduce_graph" in line:
            parallel = False
        if parallel:
            if line.endswith("vertices"):
                block_num = int(words[0][:-1])
                block_sizes[block_num] = int(words[1])
            if "Total time spent applying reductions  : " in line:
                parallel_time[num_blocks] += float(words[6])
            if "Kernel size after parallel run: " in line:
                parallel_size[num_blocks] += int(words[5])
            if "Number of isolated clique reductions" in line:
                block_num = int(words[0][:-1])
                current_block_num_isolated_cliques = int(words[-1])
                num_isolated_clique_reductions[num_blocks][block_num] += current_block_num_isolated_cliques
            if "Number of vertex fold reductions:" in line:
                block_num = int(words[0][:-1])
                current_block_num_vertex_folds = int(words[-1])
                num_vertex_fold_reductions[num_blocks][block_num] += current_block_num_vertex_folds
            if "Number of twin reductions (removed)" in line:
                block_num = int(words[0][:-1])
                current_block_num_twins_removed = int(words[-1])
                num_twins_removed[num_blocks][block_num] += current_block_num_twins_removed
            if "Number of twin reductions (folded)" in line:
                block_num = int(words[0][:-1])
                current_block_num_twins_folded = int(words[-1])
                num_twins_folded[num_blocks][block_num] += current_block_num_twins_folded
            if "Number of unconfined vertices removed" in line:
                block_num = int(words[0][:-1])
                current_block_num_unconfined_reductions = int(words[-1])
                num_unconfined_reductions[num_blocks][block_num] += current_block_num_unconfined_reductions
            if "Number of diamond reductions:" in line:
                block_num = int(words[0][:-1])
                current_block_num_diamond_reductions = int(words[-1])
                num_diamond_reductions[num_blocks][block_num] += current_block_num_diamond_reductions
            if "Number of vertices removed by LP reduction:" in line:
                num_lp_removed[num_blocks] += int(words[-1])
            if "Time spent applying reductions" in line:
                block_num = int(words[0][:-1])
                time_per_block[num_blocks][block_num] += float(words[-1])
            if "Before call to sequential reduce_graph" in line:
                parallel = False
        else:
            if "Total time spent applying reductions  : " in line:
                sequential_time[num_blocks] += float(words[6])
    sizes.sort()
    partitioning_time = getPartitionTime(inputFile, sizes[-1])
    is64 = True
    if not sizes[-1] == 64:
        is64 = False
    return (parallel_time[sizes[-1]], parallel_size[sizes[-1]], partitioning_time, is64)

ourResultsDir = sys.argv[1]
akibaResultsDir = sys.argv[2]

latexFileName = os.path.join(akibaResultsDir, "parallelComparison.tex")

latexFile = open(latexFileName, "w")

latexFile.write("\\documentclass[varwidth=\\maxdimen, border=10pt]{standalone} \n")
latexFile.write("\\usepackage[autolanguage]{numprint} \n")
latexFile.write("\\begin{document} \n")

latexFile.write("\\begin{tabular}{l|rr|rrrr|r}\n")
latexFile.write(" & \multicolumn{2}{c|}{Akiba and Iwata \cite{akiba2016branch}} & \multicolumn{4}{c|}{Our algorithm} & \\\\ \n")
latexFile.write("graph & size & kern. [s] & size & part. [s] & kern. [s] & total [s] & speedup \\\\ \n")
latexFile.write("\\hline \n")

for filename in os.listdir(ourResultsDir):
    if filename.endswith(".graph"):
        ourtime, oursize, partitioning_time, is64 = getOurTimeAndSize(os.path.join(ourResultsDir, filename))
        totalTime = ourtime + partitioning_time
        akibaFinaltime, akibaFinalsize  = getAkibaFinalTimeAndSize(akibaResultsDir, filename)
        akibaTimeToOurSize = getAkibaRuntimeForSize(akibaResultsDir, filename, oursize)

        filenameToPrint = filename.replace("-sorted", "").replace(".graph", "")
        if not is64:
            filenameToPrint += "*"
        latexFile.write(filenameToPrint + " & \\numprint{" + str(int(round(akibaFinalsize))) + "} & \\numprint{" + str(int(round(akibaFinaltime))) + "} & \\numprint{" + str(int(round(oursize))) + "} & \\numprint{" + str(int(round(partitioning_time))) + "} & \\numprint{" + str(int(round(ourtime))) + "} & \\numprint{" + str(int(round(totalTime))) + "} & \\numprint{" + str(round(akibaFinaltime / totalTime, 1)) + "} \\\\ \n" )

latexFile.write("\\end{tabular} \n")

latexFile.write("\n\n")

latexFile.write("\\begin{tabular}{l|r|r|rrr|r}\n")
latexFile.write(" &  & \multicolumn{1}{|c|}{\cite{akiba2016branch}} & \multicolumn{3}{c|}{Our algorithm} & \\\\ \n")
latexFile.write("graph & size & kern. [s] & part. [s] & kern. [s] & total [s] & speedup \\\\ \n")
latexFile.write("\\hline \n")

for filename in os.listdir(ourResultsDir):
    if filename.endswith(".graph"):
        ourtime, oursize, partitioning_time, is64 = getOurTimeAndSize(os.path.join(ourResultsDir, filename))
        totalTime = ourtime + partitioning_time
        akibaFinaltime, akibaFinalsize  = getAkibaFinalTimeAndSize(akibaResultsDir, filename)
        akibaTimeToOurSize = getAkibaRuntimeForSize(akibaResultsDir, filename, oursize)

        filenameToPrint = filename.replace("-sorted", "").replace(".graph", "")
        if not is64:
            filenameToPrint += "*"
        latexFile.write(filenameToPrint + " & \\numprint{" + str(int(round(oursize))) + "} & \\numprint{" + str(int(round(akibaTimeToOurSize))) + "} & \\numprint{" + str(int(round(partitioning_time))) + "} & \\numprint{" + str(int(round(ourtime))) + "} & \\numprint{" + str(int(round(totalTime))) + "} & \\numprint{" + str(round(akibaTimeToOurSize / totalTime, 1)) + "} \\\\ \n" )

latexFile.write("\\end{tabular} \n")
latexFile.write("\\end{document} \n")
latexFile.close()

os.system("pdflatex -output-directory " + os.path.dirname(latexFileName) + " " + latexFileName)
