import sys
import os

directory = sys.argv[1]
filename = sys.argv[2]

file = open(filename, "w")
file.write("\\documentclass[varwidth=\\maxdimen, border=10pt]{standalone} \n")
file.write("\\begin{document} \n")
file.close()

for graphDirs in os.listdir(directory):
	graphDirPath = os.path.join(directory, graphDirs)
	if os.path.isdir(graphDirPath):
		for resultsFile in os.listdir(graphDirPath):
			resultsFilePath = os.path.join(graphDirPath, resultsFile)
			if resultsFile.startswith("weight") and os.path.isfile(resultsFilePath):
				print(resultsFilePath)
				os.system("python3 draw_graphs.py " + resultsFilePath + " " + filename)

file = open(filename, "a")
file.write("\\end{document} \n")
file.close()

os.system("pdflatex -output-directory " + os.path.dirname(filename) + " " + filename)

