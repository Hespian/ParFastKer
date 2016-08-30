import sys

inputFilePath = sys.argv[1]
outputFilePath = sys.argv[2]

inputFile = open(inputFilePath, "r")
outputFile = open(outputFilePath, "w")

header = inputFile.readline()
headerWords = header.split()
numVertices = headerWords[0]
numEdges = int(headerWords[1])
outputFile.write("%%MatrixMarket matrix coordinate integer general\n")
outputFile.write(numVertices + " " + numVertices + " " + str(2 * numEdges) + "\n")
vertex = 1
numNonZeros = 0
for line in inputFile:
    neighbors = line.split()
    for neighbor in neighbors:
        outputFile.write(str(vertex) + " " + neighbor + " 1\n")
        numNonZeros += 1
    vertex += 1

if numNonZeros != 2 * numEdges:
    raise Exception(str(numNonZeros) + " != " + str(2 * numEdges))
inputFile.close()
outputFile.close()