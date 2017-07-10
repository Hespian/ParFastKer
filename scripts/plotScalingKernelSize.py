import sys
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib2tikz import save as tikz_save

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv') 
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def safe_div(x,y):
    if y == 0:
        return 0
    return x / y

def readGraphFile(inputFile, color):
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
		if "Total time spent undoing" in line:
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
			graphName = os.path.os.path.splitext(words[-1])[0]
			print(graphName)
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

	numVertices = int(numVerticesString)


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
		vertex_folds_relative[i] = safe_div(sum_vertex_folds[i], sum_vertex_folds[0]) * 100
		isolated_clique_relative[i] = safe_div(sum_isolated_clique[i], sum_isolated_clique[0]) * 100
		twins_removed_relative[i] = safe_div(sum_twins_removed[i], sum_twins_removed[0]) * 100
		twins_folded_relative[i] = safe_div(sum_twins_folded[i], sum_twins_folded[0]) * 100
		unconfined_relative[i] = safe_div(sum_unconfined[i], sum_unconfined[0]) * 100
		diamond_relative[i] = safe_div(sum_diamond[i], sum_diamond[0]) * 100
		lp_relative[i] = safe_div(num_lp_sorted[i], num_lp_sorted[0]) * 100


	parallel_times = [0] * len(sizes)
	parallel_times_rel = [0] * len(sizes)
	parallel_sizes_rel = [0] * len(sizes)
	parallel_vertices_removed = [0] * len(sizes)
	for i in range(len(sizes)):
		parallel_times[i] = parallel_time[sizes[i]]
		parallel_times_rel[i] = parallel_time[sizes[i]] / parallel_time[min(parallel_time)]
		parallel_sizes_rel[i] = parallel_size[sizes[i]] / parallel_size[min(parallel_size)]
		parallel_vertices_removed[i] = numVertices - parallel_size[sizes[i]]

	parallel_vertices_removed_rel = [0] * len(sizes)
	for i in range(len(sizes)):
		parallel_vertices_removed_rel[i] = parallel_vertices_removed[i] / parallel_vertices_removed[0]

	parallel_vertices_removed_rel = [x * 100 for x in parallel_vertices_removed_rel]
	print(parallel_vertices_removed_rel)
	plt.semilogx(sizes, parallel_vertices_removed_rel, color=color, label=graphName, basex=2, marker="x")






directory = sys.argv[1]
filename = os.path.join(directory, "kernelSize")

cmap = get_cmap(8)
fig, ax = plt.subplots()
index = np.arange(7)
width = 0.75
cmapIndex = 1

for resultsFile in os.listdir(directory):
	resultsFilePath = os.path.join(directory, resultsFile)
	if resultsFile.endswith(".graph") and os.path.isfile(resultsFilePath):
		readGraphFile(resultsFilePath, cmap(cmapIndex))
		cmapIndex += 1


plt.xlabel('Number of threads')
plt.ylabel('Number of vertices removed')

ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
plt.legend()

tikz_save(filename + '.tikz', figureheight = '\\figureheight', figurewidth = '\\figurewidth')
pp = PdfPages(filename + '.pdf')
plt.savefig(filename + '.pdf', format='pdf')



