import matplotlib.pyplot as plt
from collections import namedtuple
import numpy as np
import sys
import os

def makedir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

inputFile = sys.argv[1]
plotsDir = inputFile + "-plots"
makedir(plotsDir)

latexFileName = None
if len(sys.argv) >= 3:
	latexFileName = sys.argv[2]

file = open(inputFile, "r")
sizes = []
parallel = True
block_sizes = {}
time_per_block = {}
parallel_time = {}
parallel_size = {}
num_isolated_clique_reductions = {}
num_vertex_fold_reductions = {}
num_twins_removed = {}
num_twins_folded = {}
num_unconfined_reductions = {}
num_lp_removed = {}
sequential_time = {}
num_blocks = 0
num_reps = 0
current_rep = 0
graphName = ""
numEdgesString = ""
numVerticesString = ""
for line in file:
	words = line.split()
	if "After call to sequential reduce_graph" in line:
		if current_rep == num_reps - 1:
			for block_num in range(num_blocks):
				time_per_block[num_blocks][block_num] /= num_reps
				num_isolated_clique_reductions[num_blocks][block_num] /= num_reps
				num_vertex_fold_reductions[num_blocks][block_num] /= num_reps
				num_twins_removed[num_blocks][block_num] /= num_reps
				num_twins_folded[num_blocks][block_num] /= num_reps
				num_unconfined_reductions[num_blocks][block_num] /= num_reps
			num_lp_removed[num_blocks] /= num_reps
			parallel_time[num_blocks] /= num_reps
			parallel_size[num_blocks] /= num_reps
			sequential_time[num_blocks] /= num_reps
	if "Filename:" in line:
		graphName = words[-1]
	if "|-Nodes:" in line:
		numVerticesString = words[-1]
	if "|-Edges:" in line:
		numEdgesString = words[-1]
	if "Number of repititions" in line:
		num_reps = int(words[-1])
	if "Number of blocks:" in line:
		num_blocks = int(words[-1])
		sizes.append(num_blocks)
		block_sizes[num_blocks] = [0] * num_blocks
		time_per_block[num_blocks] = [0.0] * num_blocks
		num_isolated_clique_reductions[num_blocks] = [0] * num_blocks
		num_vertex_fold_reductions[num_blocks] = [0] * num_blocks
		num_twins_removed[num_blocks] = [0] * num_blocks
		num_twins_folded[num_blocks] = [0] * num_blocks
		num_unconfined_reductions[num_blocks] = [0] * num_blocks
		parallel_time[num_blocks] = 0
		parallel_size[num_blocks] = 0
		sequential_time[num_blocks] = 0
		num_lp_removed[num_blocks] = 0
	if "New repitition:" in line:
		parallel = True
		current_rep = int(words[-1])
	if "Before call to sequential reduce_graph" in line:
		parallel = False
	if parallel:
		if line.endswith("vertices"):
			block_num = int(words[0][:-1])
			block_sizes[block_num] = int(words[1])
		if "Total time spent applying reductions  : " in line:
			parallel_time[num_blocks] += float(words[6])
		if "Kernel size after parallel run: " in line:
			parallel_size[num_blocks] += int(words[5])
		if "Number of isolated clique reductions" in line:
			block_num = int(words[0][:-1])
			current_block_num_isolated_cliques = int(words[-1])
			num_isolated_clique_reductions[num_blocks][block_num] += current_block_num_isolated_cliques
		if "Number of vertex fold reductions:" in line:
			block_num = int(words[0][:-1])
			current_block_num_vertex_folds = int(words[-1])
			num_vertex_fold_reductions[num_blocks][block_num] += current_block_num_vertex_folds
		if "Number of twin reductions (removed)" in line:
			block_num = int(words[0][:-1])
			current_block_num_twins_removed = int(words[-1])
			num_twins_removed[num_blocks][block_num] += current_block_num_twins_removed
		if "Number of twin reductions (folded)" in line:
			block_num = int(words[0][:-1])
			current_block_num_twins_folded = int(words[-1])
			num_twins_folded[num_blocks][block_num] += current_block_num_twins_folded
		if "Number of unconfined vertices removed" in line:
			block_num = int(words[0][:-1])
			current_block_num_unconfined_reductions = int(words[-1])
			num_unconfined_reductions[num_blocks][block_num] += current_block_num_unconfined_reductions
		if "Number of vertices removed by LP reduction:" in line:
			num_lp_removed[num_blocks] += int(words[-1])
		if "Time spent applying reductions" in line:
			block_num = int(words[0][:-1])
			time_per_block[num_blocks][block_num] += float(words[-1])
		if "Before call to sequential reduce_graph" in line:
			parallel = False
	else:
		if "Total time spent applying reductions  : " in line:
			sequential_time[num_blocks] += float(words[6])

file.close()

sizes.sort()
print("Sizes:")
print(sizes)
print("Repititions:")
print(num_reps)


sum_vertex_folds = [0] * len(sizes)
sum_isolated_clique = [0] * len(sizes)
sum_twins_removed = [0] * len(sizes)
sum_twins_folded = [0] * len(sizes)
sum_unconfined = [0] * len(sizes)
num_lp_sorted = [0] * len(sizes)
for i in range(len(sizes)):
	sum_vertex_folds[i] = sum(num_vertex_fold_reductions[sizes[i]])
	sum_isolated_clique[i] = sum(num_isolated_clique_reductions[sizes[i]])
	sum_twins_removed[i] = sum(num_twins_removed[sizes[i]])
	sum_twins_folded[i] = sum(num_twins_folded[sizes[i]])
	sum_unconfined[i] = sum(num_unconfined_reductions[sizes[i]])
	num_lp_sorted[i] = num_lp_removed[sizes[i]]
index = np.arange(len(sizes))
width = 0.75
plt.figure(1, figsize=(10, 10), dpi=500)
plt.subplot(611)
plt.bar(index, sum_vertex_folds, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number vertices removed by vertex fold reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 0.7)
plt.subplot(612)
plt.bar(index, sum_isolated_clique, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number of vertices removed isolated clique reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 0.7)
plt.subplot(613)
plt.bar(index, sum_twins_removed, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number vertices removed of twin reduction (removed)')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 0.7)
plt.subplot(614)
plt.bar(index, sum_twins_folded, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number of vertices removed by twin reduction (folded)')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 0.7)
plt.subplot(615)
plt.bar(index, sum_unconfined, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number of vertices removed by unconfined reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 0.7)
plt.subplot(616)
plt.bar(index, num_lp_sorted, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number of vertices removed by lp reduction')
plt.xticks(index, sizes)
plt.savefig(os.path.join(plotsDir, "num_reductions"))
plt.clf()


parallel_times = [0] * len(sizes)
sequential_times = [0] * len(sizes)
for i in range(len(sizes)):
	parallel_times[i] = parallel_time[sizes[i]]
	sequential_times[i] = sequential_time[sizes[i]]
index = np.arange(len(sizes))
width = 0.75
plt.figure()
plt.bar(index, parallel_times, width, color="red", label='parallel')
plt.bar(index, sequential_times, width, color="#F26F6F", bottom=parallel_times, label='sequential')
plt.xlabel('# Partitions')
plt.ylabel('Run time (s)')
plt.title('Run time')
plt.xticks(index, sizes)
plt.savefig(os.path.join(plotsDir, "runtime"), bbox_inches='tight')
plt.clf()


parallel_sizes = [0] * len(sizes)
for i in range(len(sizes)):
	parallel_sizes[i] = parallel_size[sizes[i]]
index = np.arange(len(sizes))
width = 0.75
plt.figure()
plt.bar(index, parallel_sizes, width, color="red")
plt.xlabel('# Partitions')
plt.ylabel('Kernel size')
plt.title('Kernel size after parallel reductions')
plt.xticks(index, sizes)
plt.savefig(os.path.join(plotsDir, "kernel_size"))
plt.clf()


data =  [None] * len(sizes)
for i in range(len(sizes)):
	data[i] = time_per_block[sizes[i]]
plt.figure()
plt.boxplot(data, labels=sizes, showmeans=True)
plt.xlabel('# Partitions')
plt.ylabel('Run time (s)')
plt.title('Run time per partition box plots')
plt.savefig(os.path.join(plotsDir, "runtime_boxplots"))
plt.clf()

if not latexFileName == None:
	latexFile = open(latexFileName, "a")
	weight_string = os.path.basename(inputFile)
	latexFile.write(graphName + " (" + numVerticesString + ", " + numEdgesString + ") " + weight_string.replace("_","\\_") + "\n")
	latexFile.write("\\newline \n")
	latexFile.write("\\begin{tabular}{l|rrrrrrrr|r}\n")
	latexFile.write("number of blocks & parallel run time & parallel kernel size & vertex folds & isolated clique reductions & twin reductions (removed) & twin reductions (folded) & unconfined reductions & lp reductions & sequential run time \\\\ \n")
	latexFile.write("\\hline \n")
	for i in range(len(sizes)):
		latexFile.write(str(sizes[i]) + " & " + str(round(parallel_times[i], 1)) + " & " + str(int(parallel_sizes[i])) + " & " + str(int(sum_vertex_folds[i])) + " & " + str(int(sum_isolated_clique[i])) + " & " + str(int(sum_twins_removed[i])) + " & " + str(int(sum_twins_folded[i])) + " & " + str(int(sum_unconfined[i])) + " & " + str(int(num_lp_sorted[i])) + " & " + str(round(sequential_times[i], 1)) + "\\\\ \n")

	latexFile.write("\\end{tabular} \n")
	latexFile.write("\\vspace{1cm} \n\\newline \n")
	latexFile.close()
