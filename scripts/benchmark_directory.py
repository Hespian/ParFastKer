import sys
import os
import subprocess

graphDir = sys.argv[1]
numreps = sys.argv[2]

def makedir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

print("#!/bin/bash")

makeCall = "cd ../build; cmake .. -DCMAKE_BUILD_TYPE=Release; make -j48; cd ../scripts"
print(makeCall)

for graphFile in os.listdir(graphDir):
    if graphFile.endswith(".graph"):
        graphPath = os.path.join(graphDir, graphFile)
        partitionsDir = os.path.join(graphDir, "partitions", graphFile)
        for weightPartitionDir in os.listdir(partitionsDir):
            if weightPartitionDir.startswith("weight"):
                partitionsDirPath = os.path.join(partitionsDir, weightPartitionDir)
                allResultsdir = os.path.join(graphDir, "results")
                makedir(allResultsdir)
                resultsDirGraph = os.path.join(allResultsdir, graphFile)
                makedir(resultsDirGraph)
                resultsFile = os.path.join(resultsDirGraph, weightPartitionDir)
                benchmarkCall = "../build/benchmark " + graphPath + "  --partition_path=" + partitionsDirPath + " --console_log --num_reps=" + numreps + " &> " + resultsFile
                print("echo '" + benchmarkCall + "'")
                print(benchmarkCall)
