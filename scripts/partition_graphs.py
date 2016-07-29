import sys
import os
import subprocess

graphDir = sys.argv[1]
numPartitionSet = [1, 2, 4, 8, 16, 32, 48, 64]

def makedir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

print("#!/bin/bash")

def writeDegreeFile(inputFilePath, outputFilePath):
    if not os.path.exists(outputFilePath):
        inputFile = open(inputFilePath, "r")
        outputFile = open(outputFilePath, "w")

        header = inputFile.readline()
        headerWords = header.split()
        numVertices = headerWords[0]
        numEdges = headerWords[1]
        outputFile.write(numVertices + " " + numEdges + " 010\n")
        for line in inputFile:
            degree = len(line.split())
            outputFile.write(str(degree) + " " + line)

        inputFile.close()
        outputFile.close()

def writeCustomWeightFile(graphFilePath, weightsFilePath, outputFilePath):
    if not os.path.exists(outputFilePath):
        graphFile = open(graphFilePath, "r")
        weightsFile = open(weightsFilePath, "r")
        outputFile = open(outputFilePath, "w")

        header = graphFile.readline()
        headerWords = header.split()
        numVertices = headerWords[0]
        numEdges = headerWords[1]
        outputFile.write(numVertices + " " + numEdges + " 010\n")
        for line in graphFile:
            weight = weightsFile.readline().rstrip()
            outputFile.write(str(weight) + " " + line)

        graphFile.close()
        outputFile.close()

def partitionGraph(graphPath, targetDir):
    for numPartitions in numPartitionSet:
        targetFile = os.path.join(targetDir, str(numPartitions) + ".partition")
        if not os.path.exists(targetFile):
            outputDumpFile = os.path.join(targetDir, "partition_output")
            callString = "mpirun -n 16 ../../parallel_social_partitioning_package_weighted/deploy/parallel_label_compress " + graphPath + " --k=" + str(numPartitions) + " --preconfiguration=ultrafast --seed 1337 >> " + outputDumpFile
            print("echo '" + callString + "'")
            print(callString)
            outputFile = os.path.join(os.getcwd(), "tmppartition")
            mvCall = "mv " + outputFile + " " + targetFile
            print("echo '" + mvCall + "'")
            print(mvCall)


customWeightsDir = os.path.join(graphDir, "custom_weights")
customWeightFiles = []
for file in os.listdir(customWeightsDir):
    if file.endswith(".weights"):
        customWeightFiles.append(file)

for file in os.listdir(graphDir):
    if file.endswith(".graph"):
        graphPath = os.path.join(graphDir, file)
        partitionsDir = os.path.join(graphDir, "partitions", file)
        makedir(partitionsDir)
        

        # Weight one
        targetDirWeighOne = os.path.join(partitionsDir, "weight_one")
        makedir(targetDirWeighOne)
        partitionGraph(graphPath, targetDirWeighOne)

        # Weight degree
        targetDirWeighDegree = os.path.join(partitionsDir, "weight_degree")
        makedir(targetDirWeighDegree)
        weightedGraphsDir = os.path.join(graphDir, "weighted")
        makedir(weightedGraphsDir)
        weightedGraphDegreeFilePath = os.path.join(weightedGraphsDir, file) + "-weighted-degree.graph"

        writeDegreeFile(graphPath, weightedGraphDegreeFilePath)
        partitionGraph(weightedGraphDegreeFilePath, targetDirWeighDegree)

        # Custom weights
        if file + ".weights" in customWeightFiles:
            targetDirWeighCustom = os.path.join(partitionsDir, "weight_custom")
            makedir(targetDirWeighCustom)
            weightedGraphCustomFilePath = os.path.join(weightedGraphsDir, file) + "-weighted-custom.graph"
            weightsFilePath = os.path.join(customWeightsDir, file + ".weights")
            writeCustomWeightFile(graphPath, weightsFilePath, weightedGraphCustomFilePath)
            partitionGraph(weightedGraphCustomFilePath, targetDirWeighCustom)


