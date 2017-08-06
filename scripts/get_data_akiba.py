import os

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

def getAkibaFinalTimeAndSizeFromFile(filename):
    file = open(filename)
    for line in file:
        if "Final kernel size:" in line:
            words = line.split()
            size = int(words[-1])
            time = float(words[0][:-1]) / 1000.0
            result = dict()
            result["time"] = time
            result["size"] = size
            return result
    print(filename)
    raise

def getAkibaTimeAndSize(graphName, resultDir):
    time = 0.0
    size = 0
    numRuns = 0
    for filename in os.listdir(resultDir):
        if graphName in filename:
            result = getAkibaFinalTimeAndSizeFromFile(os.path.join(resultDir, filename))
            time += result["time"]
            numRuns += 1
            size = result["size"]
    time /= numRuns

    result = dict()
    result["time"] = time
    result["size"] = size
    return result

def getAkibaTimeForSize(graphName, resultDir, targetSize):
    time = 0.0
    numRuns = 0
    for filename in os.listdir(resultDir):
        if graphName in filename:
            result = getAkibaRuntimeForSizeFromFile(os.path.join(resultDir, filename), targetSize)
            time += result
            numRuns += 1
    time /= numRuns
    return time
