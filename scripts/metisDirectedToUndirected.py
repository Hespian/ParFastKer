import sys

inputFilePath = sys.argv[1]
outputFilePath = sys.argv[2]

inputFile = open(inputFilePath, "r")
outputFile = open(outputFilePath, "w")

header = inputFile.readline()
headerWords = header.split()
numVertices = int(headerWords[0])
numEdges = int(headerWords[1])

graph = [set() for i in range(numVertices + 1)]

vertex = 1
numNonZeros = 0
for line in inputFile:
    neighbors = line.split()
    neighbors = list(map(int, neighbors))
    for neighbor in neighbors:
        graph[vertex].add(neighbor)
        graph[neighbor].add(vertex)
    # print "Done with " + str(vertex)
    vertex += 1

# print "finished reading graph"
numEdges = 0
for edgeList in graph:
	numEdges += len(edgeList)

assert(numEdges % 2 == 0)
numEdges /= 2
outputFile.write(str(numVertices) + " " + str(numEdges))

for edgeList in graph:
	outputFile.write(' '.join(str(x) for x in edgeList) + "\n")

inputFile.close()
outputFile.close()