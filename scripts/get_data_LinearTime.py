import re
import os

def getLinearTimeTimeAndSizeFromFile(inputFile):
    file = open(inputFile, "r")
    numRuns = 0
    time = 0.0
    size = 0
    inputsize = 0
    for line in file:
        if "Process time" in line:
            time += float(line.split()[2][:-1]) / 1000000
            numRuns += 1
        if "n = " in line:
            words = re.split('\W+', line)
            inputsize = int(words[2])
        if "Degree_two_path MIS:" in line:
            words = re.split('\W+', line)
            size = int(words[6])
            if size == 0:
                size = inputsize
    if numRuns > 0:
        time /= numRuns
    else:
        time = -1
        size = -1
    result = dict()
    result["size"] = size
    result["time"] = time
    return result

def getLinearTimeTimeAndSize(graphName, resultsDir):
    for filename in os.listdir(resultsDir):
        if graphName in filename:
            return getLinearTimeTimeAndSizeFromFile(os.path.join(resultsDir, filename))
