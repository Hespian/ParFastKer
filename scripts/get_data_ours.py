import os

def getOurTimeAndSizeFromFiles(linearTimeOutputFile, partitioningOutputDir, ourOutputFile):
    file = open(ourOutputFile, "r")
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
    sequential_size = {}
    LP_time = {}
    Rest_time = {}
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
                num_lp_removed[num_blocks] /= num_reps
                parallel_time[num_blocks] /= num_reps
                parallel_size[num_blocks] /= num_reps
                sequential_time[num_blocks] /= num_reps
                sequential_size[num_blocks] /= num_reps
                LP_time[num_blocks] /= num_reps
                Rest_time[num_blocks] /= num_reps
        if "Filename:" in line:
            graphName = words[-1]
        if "|-Nodes:" in line:
            numVerticesString = words[-1]
        if "|-Edges:" in line:
            numEdgesString = words[-1]
        if "Number of repititions" in line:
            num_reps = int(words[-1])
        if "Number of blocks:" in line:
            num_blocks = int(words[-1])
            sizes.append(num_blocks)
            parallel_time[num_blocks] = 0
            parallel_size[num_blocks] = 0
            sequential_time[num_blocks] = 0
            sequential_size[num_blocks] = 0
            num_lp_removed[num_blocks] = 0
            LP_time[num_blocks] = 0.0
            Rest_time[num_blocks] = 0.0
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
            if "LP time: " in line:
                LP_time[num_blocks] += float(words[-1])
            if "Rest time: " in line:
                Rest_time[num_blocks] += float(words[-1])
            if "Before call to sequential reduce_graph" in line:
                parallel = False
        else:
            if "Total time spent applying reductions  : " in line:
                sequential_time[num_blocks] += float(words[6])
            if "Kernel size after sequential run: " in line:
                sequential_size[num_blocks] += int(words[5])

    file.close()
    sizes.sort()

    partitioning_times = dict()
    for size in sizes:
        partitioningOutputFile = os.path.join(partitioningOutputDir, str(size) + ".partition-log")
        file = open(partitioningOutputFile, "r")
        for line in file:
            if "avg partitioning time elapsed" in line:
                words = line.split()
                partitioning_times[size] = float(words[-1])
        file.close()

    file = open(linearTimeOutputFile, "r")
    num_LinearTime_runs = 0
    LinearTimeTime = 0.0
    for line in file:
        if "Process time" in line:
            LinearTimeTime += float(line.split()[2][:-1]) / 1000000
            num_LinearTime_runs += 1
    if num_LinearTime_runs > 0:
        LinearTimeTime /= num_LinearTime_runs
    file.close()

    max_blocks = sizes[-1]
    result = dict()
    result["lineartime_time"] = LinearTimeTime
    result["partitioning_time"] = partitioning_times[max_blocks]
    result["parallel_quasikernel_size"] = parallel_size[max_blocks]
    result["parallel_kernel_size"] = sequential_size[max_blocks]
    result["parallel_quasikernel_time"] = parallel_time[max_blocks]
    result["parallel_kernel_time"] = sequential_time[max_blocks]
    result["sequential_quasikernel_size"] = parallel_size[1]
    result["sequential_kernel_size"] = sequential_size[1]
    result["sequential_quasikernel_time"] = parallel_time[1]
    result["sequential_kernel_time"] = sequential_time[1]

    result["scaling_quasikernel_time"] = parallel_time
    result["scaling_partitioning_time"] = partitioning_times
    result["scaling_partitioning_time"][1] = 0
    result["scaling_total"] = dict()
    for i in sizes:
        result["scaling_total"][i] = result["scaling_quasikernel_time"][i] + result["scaling_partitioning_time"][i] + LinearTimeTime
    result["scaling_LP_time"] = LP_time
    result["scaling_Rest_time"] = Rest_time
    return result

def getFirstTimeWithSizeLessThan(timesAndSizes, targetSize):
    for (time, size) in timesAndSizes:
        if(size <= targetSize):
            return time
    return -100000000

def getOurTimeForSizeFromFiles(linearTimeOutputFile, partitioningOutputDir, ourOutputFile, targetSize):
    file = open(ourOutputFile, "r")
    sizes = []
    TimesAndSizes = dict()
    num_blocks = 0
    num_reps = 0
    current_rep = 0
    for line in file:
        words = line.split()
        if "Number of repititions" in line:
            num_reps = int(words[-1])
        if "Number of blocks:" in line:
            num_blocks = int(words[-1])
            sizes.append(num_blocks)
            TimesAndSizes[num_blocks] = [] * num_reps
        if "New repitition:" in line:
            current_rep = int(words[-1])
            TimesAndSizes[num_blocks].append([])
        if "finished iteration" in line:
            time = float(words[6])
            size = int(words[-1])
            tup = (time, size)
            TimesAndSizes[num_blocks][current_rep].append(tup)
        if "starts at" in line:
            time = float(words[4])
            size = int(words[-1])
            tup = (time, size)
            TimesAndSizes[num_blocks][current_rep].append(tup)
    file.close()


    sizes.sort()
    parallel_time = 0.0
    for elem in TimesAndSizes[sizes[-1]]:
        parallel_times_and_sizes = sorted(elem, key = lambda x: x[0])
        parallel_time += getFirstTimeWithSizeLessThan(parallel_times_and_sizes, targetSize)
    parallel_time /= len(TimesAndSizes[sizes[-1]])

    sequential_time = 0.0
    for elem in TimesAndSizes[1]:
        sequential_times_and_sizes = sorted(elem, key = lambda x: x[0])
        sequential_time += getFirstTimeWithSizeLessThan(sequential_times_and_sizes, targetSize)
    sequential_time /= len(TimesAndSizes[1])

    result = dict()
    results = getOurTimeAndSizeFromFiles(linearTimeOutputFile, partitioningOutputDir, ourOutputFile)
    result["parallel"] = parallel_time
    result["sequential"] = sequential_time
    return result

def getOurTimeAndSize(graphname, linearTimeDir, partitioningOutputDir, ourOutputDir):
    linearTimeOutputFile = ""
    for fileName in os.listdir(linearTimeDir):
        if graphname in fileName:
            linearTimeOutputFile = os.path.join(linearTimeDir, fileName)

    ourOutputFile = ""
    for fileName in os.listdir(ourOutputDir):
        if graphname in fileName:
            ourOutputFile = os.path.join(ourOutputDir, fileName)

    partitioningOutputDirForGraph = ""
    for dirname in os.listdir(partitioningOutputDir):
        if graphname in dirname:
            partitioningOutputDirForGraph = os.path.join(partitioningOutputDir, dirname, "weight_one_ultrafast")

    return getOurTimeAndSizeFromFiles(linearTimeOutputFile, partitioningOutputDirForGraph, ourOutputFile)

def getOurTimeForSize(graphname, linearTimeDir, partitioningOutputDir, ourOutputDir, targetSize):

    linearTimeOutputFile = ""
    for fileName in os.listdir(linearTimeDir):
        if graphname in fileName:
            linearTimeOutputFile = os.path.join(linearTimeDir, fileName)

    ourOutputFile = ""
    for fileName in os.listdir(ourOutputDir):
        if graphname in fileName:
            ourOutputFile = os.path.join(ourOutputDir, fileName)

    partitioningOutputDirForGraph = ""
    for dirname in os.listdir(partitioningOutputDir):
        if graphname in dirname:
            partitioningOutputDirForGraph = os.path.join(partitioningOutputDir, dirname, "weight_one_ultrafast")

    return getOurTimeForSizeFromFiles(linearTimeOutputFile, partitioningOutputDirForGraph, ourOutputFile, targetSize)
