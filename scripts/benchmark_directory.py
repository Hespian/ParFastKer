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

partitionCall = "python partition_graphs.py " + graphDir
print(partitionCall)
subprocess.call(partitionCall, shell=True)

makeCall = "cd ../build; cmake .. -DCMAKE_BUILD_TYPE=Release; make -j48; cd ../scripts"
print(makeCall)
subprocess.call(makeCall, shell=True)

for file in os.listdir(graphDir):
    if file.endswith(".graph"):
        graphPath = os.path.join(graphDir, file)
        partitionsDir = os.path.join(graphDir, "partitions", file)
        partitionsDirWeighOne = os.path.join(partitionsDir, "weight_one")
        allResultsdir = os.path.join(graphDir, "results")
        makedir(allResultsdir)
        resultsDirWeightOne = os.path.join(allResultsdir, "weight_one")
        makedir(resultsDirWeightOne)
        resultsFile = os.path.join(resultsDirWeightOne, file)
        benchmarkCall = "../build/benchmark " + graphPath + "  --partition_path=" + partitionsDirWeighOne + " --console_log --num_reps=" + numreps + " &> " + resultsFile
        print(benchmarkCall)
        subprocess.call(benchmarkCall, shell=True)
