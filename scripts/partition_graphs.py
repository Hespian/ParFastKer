import sys
import os
import subprocess

graphDir = sys.argv[1]
# numPartitionSet = [1, 2, 4, 8, 16, 32, 64]
# numPartitionSet = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
numPartitionSet = [1, 32]
preconfigurations = ["ultrafast"]

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

def writeTwoNeighborhoodFile(inputFilePath, outputFilePath):
    if not os.path.exists(outputFilePath):
        inputFile = open(inputFilePath, "r")

        number_nodes = 0
        number_edges = 0
        edges_counted = 0
        adjacency = []
        node = 0
        for line in inputFile:
            args = line.strip().split()
            if node == 0:
                number_nodes = int(args[0])
                number_edges = int(args[1])
                adjacency = [[] for _ in range(0, number_nodes + 1)]
                adjacency[0] = [str(number_nodes), str(number_edges)]
            else:
                adjacency[node] = args
                edges_counted += len(args)
            node += 1
        if edges_counted < number_edges:
            print "Found less edges than specified"
            sys.exit(0)
        inputFile.close()

        weight = [0] * (number_nodes + 1)
        for i in range(1, node):
            for j in adjacency[i]:
                neighbor = int(j)
                neighborDegree = len(adjacency[neighbor])
                weight[i] += neighborDegree

        outputFile = open(outputFilePath, "w")
        outputFile.write(str(number_nodes) + " " + str(number_edges) + " 010\n")
        for i in range(1, node):
            outputFile.write(str(weight[i]) + " " + ' '.join(adjacency[i]) + "\n")

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

def partitionGraph(graphPath, targetDir, preconfiguration):
    for numPartitions in numPartitionSet:
        # if ("uk-2007-05" in graphPath) and (numPartitions == 64):
        #     numPartitions = 32
        targetFile = os.path.join(targetDir, str(numPartitions) + ".partition")
        if not os.path.exists(targetFile):
            outputDumpFile = os.path.join(targetDir, "partition_output")
            callString = "mpirun -n " + "32" + " ../../parallel_social_partitioning_package_weighted/deploy/parallel_label_compress_reps " + graphPath + " --k=" + str(numPartitions) + " --seed 1337 --num_tries=5 --preconfiguration=" + preconfiguration + " &> " + outputDumpFile
            print("echo '" + callString + "'")
            print(callString)
            outputFile = os.path.join(os.getcwd(), "tmppartition")
            mvCall = "mv " + outputFile + " " + targetFile
            print("echo '" + mvCall + "'")
            print(mvCall)
            logfiletarget = targetFile + "-log"
            mvCallLogfile = "mv " + outputDumpFile + " " + logfiletarget
            print("echo '" + mvCallLogfile + "'")
            print(mvCallLogfile)

def partitionGraphLPA(graphPath, targetDir):
    for numPartitions in numPartitionSet:
        targetFile = os.path.join(targetDir, str(numPartitions) + ".partition")
        if not os.path.exists(targetFile):
            outputDumpFile = os.path.join(targetDir, "partition_output")
            callString = "../../KaHIPLPkway/deploy/label_propagation " + graphPath + " --k=" + str(numPartitions) +  " &> " + outputDumpFile
            print("echo '" + callString + "'")
            print(callString)
            outputFile = os.path.join(os.getcwd(), "tmpclustering")
            mvCall = "mv " + outputFile + " " + targetFile
            print("echo '" + mvCall + "'")
            print(mvCall)
            logfiletarget = targetFile + "-log"
            mvCallLogfile = "mv " + outputDumpFile + " " + logfiletarget
            print("echo '" + mvCallLogfile + "'")
            print(mvCallLogfile)

customWeightsDir = os.path.join(graphDir, "custom_weights")
customWeightFiles = []
for file in os.listdir(customWeightsDir):
    if file.endswith(".weights"):
        customWeightFiles.append(file)
'''
for file in os.listdir(graphDir):
    if file.endswith(".graph"):
        graphPath = os.path.join(graphDir, file)
        partitionsDir = os.path.join(graphDir, "partitions", file)
        makedir(partitionsDir)


        # Weight one                                                                                                                                                                                       
        targetDirWeighOne = os.path.join(partitionsDir, "weight_one_LPA")
        makedir(targetDirWeighOne)
        partitionGraphLPA(graphPath, targetDirWeighOne)

        # Weight degree                                                                                                                                                                                     
        weightedGraphsDir = os.path.join(graphDir, "weighted")
        makedir(weightedGraphsDir)
        weightedGraphDegreeFilePath = os.path.join(weightedGraphsDir, file) + "-weighted-degree.graph"

        writeDegreeFile(graphPath, weightedGraphDegreeFilePath)

        targetDirWeighDegree = os.path.join(partitionsDir, "weight_degree_LPA")
        makedir(targetDirWeighDegree)
        partitionGraphLPA(weightedGraphDegreeFilePath, targetDirWeighDegree)

'''
for file in os.listdir(graphDir):
    if file.endswith(".graph"):
        graphPath = os.path.join(graphDir, file)
        partitionsDir = os.path.join(graphDir, "partitions", file)
        makedir(partitionsDir)
        

        # Weight one
        for preconfiguration in preconfigurations:
            targetDirWeighOne = os.path.join(partitionsDir, "weight_one_" + preconfiguration)
            makedir(targetDirWeighOne)
            partitionGraph(graphPath, targetDirWeighOne, preconfiguration)
        continue
        # Weight degree
        weightedGraphsDir = os.path.join(graphDir, "weighted")
        makedir(weightedGraphsDir)
        weightedGraphDegreeFilePath = os.path.join(weightedGraphsDir, file) + "-weighted-degree.graph"

        writeDegreeFile(graphPath, weightedGraphDegreeFilePath)

        for preconfiguration in preconfigurations:
            targetDirWeighDegree = os.path.join(partitionsDir, "weight_degree_" + preconfiguration)
            makedir(targetDirWeighDegree)
            partitionGraph(weightedGraphDegreeFilePath, targetDirWeighDegree, preconfiguration)

        # Weight two neighborhood size
        # targetDirWeighTwoNeighborhood = os.path.join(partitionsDir, "weight_2_neighborhood")
        # makedir(targetDirWeighTwoNeighborhood)
        # weightedGraphTwoNeighborhoodFilePath = os.path.join(weightedGraphsDir, file) + "-weighted-2-neighborhood.graph"
        # writeTwoNeighborhoodFile(graphPath, weightedGraphTwoNeighborhoodFilePath)
        # partitionGraph(weightedGraphTwoNeighborhoodFilePath, targetDirWeighTwoNeighborhood)

        # Custom weights
        if file + ".weights" in customWeightFiles:
            targetDirWeighCustom = os.path.join(partitionsDir, "weight_custom")
            makedir(targetDirWeighCustom)
            weightedGraphCustomFilePath = os.path.join(weightedGraphsDir, file) + "-weighted-custom.graph"
            weightsFilePath = os.path.join(customWeightsDir, file + ".weights")
            writeCustomWeightFile(graphPath, weightsFilePath, weightedGraphCustomFilePath)
            partitionGraph(weightedGraphCustomFilePath, targetDirWeighCustom)


