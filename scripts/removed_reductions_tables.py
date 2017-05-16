import os
import sys
import operator

inputDir = sys.argv[1]

latexFileName = os.path.join(inputDir, "tables.tex")

latexFile = open(latexFileName, 'w')


# total_time, partitioned_time, kernel_size
def readGraphResults(graphResultFile):
	file = open(graphResultFile, "r")
	partition_type = ""
	num_reps = 0;
	total_time = 0.0
	partitioned_time = 0.0;
	graphName = ""
	current_rep = 0
	kernel_size = 0
	num_blocks = 0

	results = dict()

	for line in file:
		words = line.split()
		if "weight_" in line:
			total_time = 0.0
			partitioned_time = 0.0
			kernel_size = 0
			partition_type = words[0]
		if "Total time spent undoing" in line:
			if current_rep == num_reps - 1:
				total_time /= num_reps
				partitioned_time /= num_reps
				kernel_size /= num_reps
				if partition_type == "weight_one_ultrafast" and num_blocks == 16:
					results = (total_time, partitioned_time, round(kernel_size))
		if "Filename:" in line:
			graphName = words[-1]
		if "Number of repititions" in line:
			num_reps = int(words[-1])
		if "New repitition:" in line:
			current_rep = int(words[-1])
		if "Total time spent applying reductions  : " in line:
			total_time += float(words[6])
		if "Kernel size after parallel run: " in line:
			kernel_size += int(words[5])
		if "Rest time" in line:
			partitioned_time += float(words[-1])
		if "Number of blocks:" in line:
			num_blocks = int(words[-1])
			total_time = 0.0
			partitioned_time = 0.0;
			current_rep = 0
			kernel_size = 0
	file.close()
	return results


def addTable(graphResults):
	latexFile.write("\\begin{tabular}{l|rr}\n")
	latexFile.write("missing reduction & reduction time & kernel size \\\\ \n")
	latexFile.write("\\hline \n")

	for missingReduction, results in iter(sorted(graphResults.items())):
		latexFile.write(missingReduction.replace("_","\\_") + " & " + str(round(results[0] * 100)) + " & " + str(round(results[2] * 100)) + " \\\\ \n")

	latexFile.write("\\end{tabular} \n")
	latexFile.write("\\vspace{1cm} \n\\newline \n")

def addSummaryTable(avg_results):
	latexFile.write("\\begin{tabular}{l|rr}\n")
	latexFile.write("configuration & reduction time & kernel size \\\\ \n")
	latexFile.write("\\hline \n")

	for missingReduction, results in iter(sorted(avg_results.items())):
		latexFile.write(missingReduction.replace("_"," ") + " & " + '{0:+}'.format(round(results[0] * 100 - 100, 1)) + "\% & " + '{0:+}'.format(round(results[2] * 100 - 100, 1)) + "\% \\\\ \n")

	latexFile.write("\\end{tabular} \n")
	latexFile.write("\\vspace{1cm} \n\\newline \n")


results = dict()
for dirname in os.listdir(inputDir):
	dirResults = dict()
	dirPath = os.path.join(inputDir, dirname)
	if not os.path.isdir(dirPath):
		continue
	for filename in os.listdir(dirPath):
	    if filename.endswith(".graph"): 
	        dirResults[filename] = readGraphResults(os.path.join(dirPath, filename))
	        continue
	    else:
	        continue
	results[dirname] = dirResults


base_case = "all_reductions"
rel_results = dict()
for dirname, graphdict in results.items():
	for graphName, graphResults in graphdict.items():
		if not graphName in rel_results:
			rel_results[graphName] = dict()
		rel_per_graph = (graphResults[0] / results[base_case][graphName][0], graphResults[1] / results[base_case][graphName][1], graphResults[2] / results[base_case][graphName][2])
		rel_results[graphName][dirname] = rel_per_graph

avg_results = dict()
for graph, resultDict in rel_results.items():
	for reduction, reductionResults in resultDict.items():
		oldvalue = avg_results.get(reduction, (0,0,0))
		avg_results[reduction] = tuple(map(operator.add, oldvalue, reductionResults))

num_graphs = len(rel_results)
for reduction, resultSum in avg_results.items():
	avg_results[reduction] = (resultSum[0] / num_graphs, resultSum[1] / num_graphs, resultSum[2] / num_graphs)

latexFile.write("\\documentclass[varwidth=\\maxdimen, border=10pt]{standalone} \n")
latexFile.write("\\begin{document} \n")

for graph, graphResults in rel_results.items():
	latexFile.write(graph.replace("_","\\_") + "\n")
	latexFile.write("\\newline \n")
	addTable(graphResults)

latexFile.write("Summary relative \n")
latexFile.write("\\newline \n")
addSummaryTable(avg_results)

latexFile.write("\\end{document} \n")
latexFile.close()

os.system("pdflatex -output-directory " + os.path.dirname(latexFileName) + " " + latexFileName)


