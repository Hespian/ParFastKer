import matplotlib.pyplot as plt
from collections import namedtuple
import numpy as np
import sys
import os
import operator

def makedir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

inputFile = sys.argv[1]
plotsDir = inputFile + "-plots"
makedir(plotsDir)

latexFileName = os.path.join(plotsDir, "table.tex")


file = open(inputFile, "r")
thread_counts = []
loading_time = dict()
init_time = dict()
maximum_matching_time = dict()
vertex_marking_time = dict()
apply_results_time = dict()
graphName = ""
numEdgesString = ""
numVerticesString = ""
current_rep = 0
num_threads = 0
num_reps = 0
for line in file:
	words = line.split()
	if "Total time spent undoing  reductions" in line:
		if current_rep == num_reps - 1:
			loading_time[num_threads] = loading_time[num_threads] / num_reps 
			init_time[num_threads] = init_time[num_threads] / num_reps 
			maximum_matching_time[num_threads] = maximum_matching_time[num_threads] / num_reps 
			vertex_marking_time[num_threads] = vertex_marking_time[num_threads] / num_reps 
			apply_results_time[num_threads] = apply_results_time[num_threads] / num_reps 
	if "Filename:" in line:
		graphName = words[-1]
	if "|-Nodes:" in line:
		numVerticesString = words[-1]
	if "|-Edges:" in line:
		numEdgesString = words[-1]
	if "Number of repititions" in line:
		num_reps = int(words[-1])
	if "num threads:" in line:
		num_threads = int(words[-1])
		if num_threads not in thread_counts:
			thread_counts.append(num_threads)
	if "New repitition:" in line:
		current_rep = int(words[-1])
	if "Time for loading the graph:" in line:
		time = float(words[-1])
		loading_time[num_threads] = loading_time.get(num_threads, 0.0) + time
	if "Time for KarpSipserInit:" in line:
		time = float(words[-1])
		init_time[num_threads] = init_time.get(num_threads, 0.0) + time
	if "Time for MS_BFS_Graft:" in line:
		time = float(words[-1])
		maximum_matching_time[num_threads] = maximum_matching_time.get(num_threads, 0.0) + time
	if "Time for MarkReachableVertices:" in line:
		time = float(words[-1])
		vertex_marking_time[num_threads] = vertex_marking_time.get(num_threads, 0.0) + time
	if "Time for applying result:" in line:
		time = float(words[-1])
		apply_results_time[num_threads] = apply_results_time.get(num_threads, 0.0) + time

file.close()

thread_counts.sort()

loading_time_sorted = [0.0] * len(thread_counts)
init_time_sorted = [0.0] * len(thread_counts)
maximum_matching_time_sorted = [0.0] * len(thread_counts)
vertex_marking_time_sorted = [0.0] * len(thread_counts)
apply_results_time_sorted = [0.0] * len(thread_counts)
for i in range(len(thread_counts)):
	loading_time_sorted[i] = loading_time[thread_counts[i]]
	init_time_sorted[i] = init_time[thread_counts[i]]
	maximum_matching_time_sorted[i] = maximum_matching_time[thread_counts[i]]
	vertex_marking_time_sorted[i] = vertex_marking_time[thread_counts[i]]
	apply_results_time_sorted[i] = apply_results_time[thread_counts[i]]
index = np.arange(len(thread_counts))
width = 0.75
# plt.figure(1, figsize=(10, 10), dpi=500)
plt.bar(index, loading_time_sorted, width, color="red", label='Loading graph')
plt.bar(index, init_time_sorted, width, bottom=loading_time_sorted, color="green", label='Karp Sipser init')
sum_times_1 = list(map(operator.add, loading_time_sorted, init_time_sorted))
plt.bar(index, maximum_matching_time_sorted, width, bottom=sum_times_1, color="blue", label='Maximum matching')
sum_times_2 = list(map(operator.add, sum_times_1, maximum_matching_time_sorted))
plt.bar(index, vertex_marking_time_sorted, width, bottom=sum_times_2, color="yellow", label='Mark reachable vertices')
sum_times_3 = list(map(operator.add, sum_times_2, vertex_marking_time_sorted))
plt.bar(index, apply_results_time_sorted, width, bottom=sum_times_3, color="orange", label='Apply results')
plt.xlabel('# Threads')
plt.ylabel('Time [s]')
plt.title('Time for applying the LP reduction in parallel')
plt.xticks(index, thread_counts)
plt.legend()
plt.savefig(os.path.join(plotsDir, "LP_reduction_times"))
plt.clf()




if not latexFileName == None:
	latexFile = open(latexFileName, "w")
	latexFile.write("\\documentclass[varwidth=\\maxdimen, border=10pt]{standalone} \n")
	latexFile.write("\\begin{document} \n")
	latexFile.write("\\begin{tabular}{l|rrrrr}\n")
	latexFile.write("Number of threads & Loading graph & Karp Sipser init & Maximum matching & Mark reachable vertices & Apply results \\\\ \n")
	latexFile.write("\\hline \n")
	latexFile.write(str(thread_counts[0]) + "[seconds] & " + str(round(loading_time_sorted[0], 2)) + " & " + str(round(init_time_sorted[0], 2)) + " & " + str(round(maximum_matching_time_sorted[0], 2)) + " & " + str(round(vertex_marking_time_sorted[0], 2)) + " & " + str(round(apply_results_time_sorted[0], 2)) + "\\\\ \n")
	latexFile.write("\\hline \n")
	for i in range(1, len(thread_counts)):
		latexFile.write(str(thread_counts[i]) + "[speedup] & " + str(round(loading_time_sorted[0] / loading_time_sorted[i], 2)) + " & " + str(round(init_time_sorted[0] / init_time_sorted[i], 2)) + " & " + str(round(maximum_matching_time_sorted[0] / maximum_matching_time_sorted[i], 2)) + " & " + str(round(vertex_marking_time_sorted[0] / vertex_marking_time_sorted[i], 2)) + " & " + str(round(apply_results_time_sorted[0] / apply_results_time_sorted[i], 2)) + "\\\\ \n")

	latexFile.write("\\end{tabular} \n")
	latexFile.write("\\vspace{1cm} \n")
	latexFile.write("\\end{document} \n")
	latexFile.close()
	os.system("pdflatex -output-directory " + os.path.dirname(latexFileName) + " " + latexFileName)
