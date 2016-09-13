import matplotlib.pyplot as plt
from collections import namedtuple
import numpy as np
import sys
import os
import operator
import matplotlib.cm as cmx
import matplotlib.colors as colors

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

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
num_diamond_reductions = {}
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
				num_diamond_reductions[num_blocks][block_num] /= num_reps
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
		num_diamond_reductions[num_blocks] = [0] * num_blocks
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
		if "Number of diamond reductions:" in line:
			block_num = int(words[0][:-1])
			current_block_num_diamond_reductions = int(words[-1])
			num_diamond_reductions[num_blocks][block_num] += current_block_num_diamond_reductions
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
sum_diamond = [0] * len(sizes)
num_lp_sorted = [0] * len(sizes)
for i in range(len(sizes)):
	sum_vertex_folds[i] = sum(num_vertex_fold_reductions[sizes[i]])
	sum_isolated_clique[i] = sum(num_isolated_clique_reductions[sizes[i]])
	sum_twins_removed[i] = sum(num_twins_removed[sizes[i]])
	sum_twins_folded[i] = sum(num_twins_folded[sizes[i]])
	sum_unconfined[i] = sum(num_unconfined_reductions[sizes[i]])
	sum_diamond[i] = sum(num_diamond_reductions[sizes[i]])
	num_lp_sorted[i] = num_lp_removed[sizes[i]]

vertex_folds_relative = [0] * len(sizes)
isolated_clique_relative = [0] * len(sizes)
twins_removed_relative = [0] * len(sizes)
twins_folded_relative= [0] * len(sizes)
unconfined_relative = [0] * len(sizes)
diamond_relative = [0] * len(sizes)
lp_relative = [0] * len(sizes)
for i in range(len(sizes)):
	vertex_folds_relative[i] = (sum_vertex_folds[i] / sum_vertex_folds[0]) * 100
	isolated_clique_relative[i] = (sum_isolated_clique[i] / sum_isolated_clique[0]) * 100
	twins_removed_relative[i] = (sum_twins_removed[i] / sum_twins_removed[0]) * 100
	twins_folded_relative[i] = (sum_twins_folded[i] / sum_twins_folded[0]) * 100
	unconfined_relative[i] = (sum_unconfined[i] / sum_unconfined[0]) * 100
	diamond_relative[i] = (sum_diamond[i] / sum_diamond[0]) * 100
	lp_relative[i] = (num_lp_sorted[i] / num_lp_sorted[0]) * 100

index = np.arange(len(sizes))
width = 0.75
plt.figure(1, figsize=(10, 10), dpi=500)
plt.subplot(711)
plt.bar(index, vertex_folds_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed by vertex fold reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 1.3)
plt.subplot(712)
plt.bar(index, isolated_clique_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed isolated clique reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 1.3)
plt.subplot(713)
plt.bar(index, twins_removed_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed of twin reduction (removed)')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 1.3)
plt.subplot(714)
plt.bar(index, twins_folded_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed by twin reduction (folded)')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 1.3)
plt.subplot(715)
plt.bar(index, unconfined_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed by unconfined reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 1.3)
plt.subplot(716)
plt.bar(index, lp_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed by lp reduction')
plt.xticks(index, sizes)
plt.subplots_adjust(hspace = 1.3)
plt.subplot(717)
plt.bar(index, diamond_relative, width, color="red")
plt.xlabel('# Threads')
plt.ylabel('% Vertices')
plt.title('Vertices removed by diamond reduction')
plt.xticks(index, sizes)
plt.savefig(os.path.join(plotsDir, "num_reductions_seperate"))
plt.clf()

cmap = get_cmap(7)
plt.figure()
index = np.arange(len(sizes))
width = 0.75
plt.bar(index, sum_vertex_folds, width, color=cmap(0), label="vertex fold")
plt.bar(index, sum_isolated_clique, width, bottom=sum_vertex_folds, color=cmap(1), label="isolated clique")
sum_reductions_1 = list(map(operator.add, sum_vertex_folds, sum_isolated_clique))
plt.bar(index, sum_twins_removed, width, bottom=sum_reductions_1, color=cmap(2), label="twin (removed)")
sum_reductions_2 = list(map(operator.add, sum_reductions_1, sum_twins_removed))
plt.bar(index, sum_twins_folded, width, bottom=sum_reductions_2, color=cmap(3), label="twin (folded)")
sum_reductions_3 = list(map(operator.add, sum_reductions_2, sum_twins_folded))
plt.bar(index, sum_unconfined, width, bottom=sum_reductions_3, color=cmap(4), label="unconfined")
sum_reductions_4 = list(map(operator.add, sum_reductions_3, sum_unconfined))
plt.bar(index, sum_diamond, width, bottom=sum_reductions_4, color=cmap(5), label="diamond")
sum_reductions_5 = list(map(operator.add, sum_reductions_4, sum_diamond))
plt.bar(index, num_lp_sorted, width, bottom=sum_reductions_5, color=cmap(6), label="lp")
sum_reductions_6 = list(map(operator.add, sum_reductions_5, num_lp_sorted))
plt.xlabel('# Threads')
plt.ylabel('# Vertices')
plt.title('Number of vertices removed by reductions')
plt.xticks(index, sizes)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(os.path.join(plotsDir, "num_reductions"), bbox_inches='tight')
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
plt.xlabel('# Threads')
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
	latexFile.write("\\begin{tabular}{l|rrrrrrrrr|r}\n")
	latexFile.write("number of blocks & parallel run time & parallel kernel size & vertex folds & isolated clique reductions & twin reductions (removed) & twin reductions (folded) & unconfined reductions & lp reductions & diamond reductions & sequential run time \\\\ \n")
	latexFile.write("\\hline \n")
	latexFile.write(str(sizes[0]) + " & " + str(round(parallel_times[0], 1)) + " & " + str(int(parallel_sizes[0])) + " & " + str(int(sum_vertex_folds[0])) + " & " + str(int(sum_isolated_clique[0])) + " & " + str(int(sum_twins_removed[0])) + " & " + str(int(sum_twins_folded[0])) + " & " + str(int(sum_unconfined[0])) + " & " + str(int(num_lp_sorted[0])) + " & " + str(int(sum_diamond[0])) + " & " + str(round(sequential_times[0], 1)) + "\\\\ \n")
	for i in range(1, len(sizes)):
		latexFile.write(str(sizes[i]) + "(\\%) & " + str(round((parallel_times[i] / parallel_times[0]) * 100, 1)) + " & " + str(int((parallel_sizes[i] / parallel_sizes[0]) * 100)) + " & " + str(int(vertex_folds_relative[i])) + " & " + str(int(isolated_clique_relative[i])) + " & " + str(int(twins_removed_relative[i])) + " & " + str(int(twins_folded_relative[i])) + " & " + str(int(unconfined_relative[i])) + " & " + str(int(lp_relative[i])) + " & " + str(int(diamond_relative[i])) + " & " + str(round((sequential_times[i] / sequential_times[0]) * 100, 1)) + "\\\\ \n")

	latexFile.write("\\end{tabular} \n")
	latexFile.write("\\vspace{1cm} \n\\newline \n")
	latexFile.close()
