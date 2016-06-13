import matplotlib.pyplot as plt
from collections import namedtuple
import numpy as np
import pickle

base_path = "../../results/uk-2002.graph"
partition_sizes = [1, 2, 4, 8, 16, 32, 48, 64]
modes = ["kahip", "parallel_kahip", "lpa"]

GraphData = namedtuple('GraphData', 'preprocessing_time fold_percentage parallel_time parallel_size sequential_time sequential_size num_isolated_clique_reductions num_vertex_fold_reductions isolated_clique_start_times vertex_fold_start_times')

def read_file(file):
	file = open(file, "r")
	fold_percentage = 0
	preprocessing_time = 0.0
	for line in file:
		words = line.split()
		if "|-Nodes:" in line:
			graph_size = int(words[-1])
		if "Sequential took" in line:
			sequential_time = float(words[2])
		if "Parallel took" in line:
			parallel_time = float(words[2])
		if "Preprocessing took" in line:
			preprocessing_time = float(words[2])
		if "nodes are considered for vertex folding" in line:
			fold_percentage = float(words[9][1:-2])
		if "Kernel size after parallel run:" in line:
			parallel_size = int(words[5])
		if "Kernel size:" in line:
			sequential_size = int(words[-1])
		if "Isolated clique times:" in line:
			words = line[22:-3].split(',')
			isolated_clique_times = list(map(float, words))
		if "Vertex fold times:" in line:
			words = line[18:-3].split(',')
			vertex_fold_times = list(map(float, words))
		if "Isolated clique start times:" in line:
			words = line[28:-3].split(',')
			isolated_clique_start_times = list(map(float, words))
		if "Vertex fold start times:" in line:
			words = line[24:-3].split(',')
			vertex_fold_start_times = list(map(float, words))
	file.close()

	num_isolated_clique_reductions = [0] * (len(isolated_clique_start_times) - 1)
	for i in range(0, len(isolated_clique_start_times) - 1):
		lower_time = isolated_clique_start_times[i]
		upper_time = isolated_clique_start_times[i + 1]
		times_in_range = [x for x in isolated_clique_times if x >= lower_time and x < upper_time]
		num_isolated_clique_reductions[i] = len(times_in_range)

	assert sum(num_isolated_clique_reductions) == len(isolated_clique_times), "Counting went wrong"

	num_vertex_fold_reductions = [0] * (len(vertex_fold_start_times) - 1)
	for i in range(0, len(vertex_fold_start_times) - 1):
		lower_time = vertex_fold_start_times[i]
		upper_time = vertex_fold_start_times[i + 1]
		times_in_range = [x for x in vertex_fold_times if x >= lower_time and x < upper_time]
		num_vertex_fold_reductions[i] = len(times_in_range)

	assert sum(num_vertex_fold_reductions) == len(vertex_fold_times), "Counting went wrong"

	# print("Graph size: 				" + str(graph_size))
	# print("Preprocessing time: 			" + str(preprocessing_time))
	# print("Fold percentage: 			" + str(fold_percentage))
	# print("Parallel time: 				" + str(parallel_time))
	# print("Parallel size: 				" + str(parallel_size))
	# print("Isolated clique reductions:		" + str(len(isolated_clique_times)))
	# print("Isolated clique tries:			" + str(len(isolated_clique_start_times)))
	# print("Vertex fold reductions: 		" + str(len(vertex_fold_times)))
	# print("Vertex fold tries: 			" + str(len(vertex_fold_start_times)))
	# print("Sequential time: 			" + str(sequential_time))
	# print("Sequential size: 			" + str(sequential_size))

	return GraphData(preprocessing_time, fold_percentage, parallel_time, parallel_size, sequential_time, sequential_size, num_isolated_clique_reductions, num_vertex_fold_reductions, isolated_clique_start_times, vertex_fold_start_times)

def draw_reductions_over_time(num_partitions, isolated_clique_start_times, num_isolated_clique_reductions, vertex_fold_start_times, num_vertex_fold_reductions):
	plt.yscale('log')
	plt.xlabel('Time')
	plt.ylabel('Number of reductions')
	plt.title('Reductions over time eduction')
	plt.plot(isolated_clique_start_times[:-1], num_isolated_clique_reductions, 'rx', label='Isolated clique reduction')
	plt.plot(vertex_fold_start_times[:-1], num_vertex_fold_reductions, 'bx', label='Vertex fold reduction')
	legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
	plt.savefig("graphics/reductions_over_time_" + str(num_partitions), bbox_extra_artists=(legend,), bbox_inches='tight')
	plt.clf()

# graphs = {}
# for num_partitions in partition_sizes:
# 	for mode in modes:
# 		graphs[(num_partitions, mode)] = read_file(base_path + "-" + str(num_partitions) + "-" + mode)

# pickle.dump( graphs, open( "graphics/graphs-no-preprocessing.p", "wb" ) )
# exit(0)
graphs = pickle.load( open( "graphics/graphs-no-preprocessing.p", "rb" ) )

for num_partitions in partition_sizes:
	draw_reductions_over_time(num_partitions, graphs[(num_partitions, "parallel_kahip")].isolated_clique_start_times, graphs[(num_partitions, "parallel_kahip")].num_isolated_clique_reductions, graphs[(num_partitions, "parallel_kahip")].vertex_fold_start_times, graphs[(num_partitions, "parallel_kahip")].num_vertex_fold_reductions)

# fold_percentages_kahip = [0] * len(partition_sizes)
# fold_percentages_parallel_kahip = [0] * len(partition_sizes)
# fold_percentages_lpa = [0] * len(partition_sizes)
# for i in range(len(partition_sizes)):
# 	fold_percentages_kahip[i] = graphs[(partition_sizes[i], "kahip")].fold_percentage
# 	fold_percentages_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].fold_percentage
# 	fold_percentages_lpa[i] = graphs[(partition_sizes[i], "lpa")].fold_percentage
# index = np.arange(len(partition_sizes))
# width = 0.25
# plt.bar(index, fold_percentages_kahip, width, color="blue", label='kahip')
# plt.bar(index + width, fold_percentages_parallel_kahip, width, color="red", label='parallel kahip')
# plt.bar(index + 2 * width, fold_percentages_lpa, width, color="green", label='lpa')
# plt.xlabel('Partitions')
# plt.ylabel('Vertices considered for fold reduction (%)')
# plt.title('Vertices considered for fold reduction')
# plt.xticks(index + 1.5 * width, partition_sizes)
# plt.legend(loc='upper right')
# plt.savefig("graphics/fold_percentage")
# plt.clf()

sum_vertex_fold_kahip = [0] * len(partition_sizes)
sum_isolated_clique_kahip = [0] * len(partition_sizes)
sum_vertex_fold_parallel_kahip = [0] * len(partition_sizes)
sum_isolated_clique_parallel_kahip = [0] * len(partition_sizes)
sum_vertex_fold_lpa = [0] * len(partition_sizes)
sum_isolated_clique_lpa = [0] * len(partition_sizes)
for i in range(len(partition_sizes)):
	sum_vertex_fold_kahip[i] = sum(graphs[(partition_sizes[i], "kahip")].num_vertex_fold_reductions)
	sum_vertex_fold_parallel_kahip[i] = sum(graphs[(partition_sizes[i], "parallel_kahip")].num_vertex_fold_reductions)
	sum_vertex_fold_lpa[i] = sum(graphs[(partition_sizes[i], "lpa")].num_vertex_fold_reductions)
	sum_isolated_clique_kahip[i] = sum(graphs[(partition_sizes[i], "kahip")].num_isolated_clique_reductions)
	sum_isolated_clique_parallel_kahip[i] = sum(graphs[(partition_sizes[i], "parallel_kahip")].num_isolated_clique_reductions)
	sum_isolated_clique_lpa[i] = sum(graphs[(partition_sizes[i], "lpa")].num_isolated_clique_reductions)
index = np.arange(len(partition_sizes))
width = 0.25
plt.figure(1)
plt.subplot(211)
plt.bar(index, sum_vertex_fold_kahip, width, color="blue", label='kahip')
plt.bar(index + width, sum_vertex_fold_parallel_kahip, width, color="red", label='parallel kahip')
plt.bar(index + 2 * width, sum_vertex_fold_lpa, width, color="green", label='lpa')
plt.xlabel('# Partitions')
plt.ylabel('# Reductions')
plt.title('Number of vertex fold reductions performed in parallel')
plt.xticks(index + 1.5 * width, partition_sizes)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.06), fancybox=True, shadow=True, ncol=3)
plt.subplots_adjust(hspace = 0.5)
plt.subplot(212)
plt.bar(index, sum_isolated_clique_kahip, width, color="blue", label='kahip')
plt.bar(index + width, sum_isolated_clique_parallel_kahip, width, color="red", label='parallel kahip')
plt.bar(index + 2 * width, sum_isolated_clique_lpa, width, color="green", label='lpa')
plt.xlabel('# Partitions')
plt.ylabel('# Reductions')
plt.title('Number of isolated clique reductions performed in parallel')
plt.xticks(index + 1.5 * width, partition_sizes)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.06), fancybox=True, shadow=True, ncol=3)
plt.savefig("graphics/num_reductions_parallel")
plt.clf()


parallel_time_kahip = [0] * len(partition_sizes)
sequential_time_kahip = [0] * len(partition_sizes)
parallel_time_parallel_kahip = [0] * len(partition_sizes)
sequential_time_parallel_kahip = [0] * len(partition_sizes)
parallel_time_lpa = [0] * len(partition_sizes)
sequential_time_lpa = [0] * len(partition_sizes)
for i in range(len(partition_sizes)):
	parallel_time_kahip[i] = graphs[(partition_sizes[i], "kahip")].parallel_time
	parallel_time_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].parallel_time
	parallel_time_lpa[i] = graphs[(partition_sizes[i], "lpa")].parallel_time
	sequential_time_kahip[i] = graphs[(partition_sizes[i], "kahip")].sequential_time
	sequential_time_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].sequential_time
	sequential_time_lpa[i] = graphs[(partition_sizes[i], "lpa")].sequential_time
index = np.arange(len(partition_sizes))
width = 0.25
plt.bar(index, parallel_time_kahip, width, color="blue", label='parallel: kahip')
plt.bar(index, sequential_time_kahip, width, color="#6ECBEB", bottom=parallel_time_kahip, label='sequential: kahip')
plt.bar(index + width, parallel_time_parallel_kahip, width, color="red", label='parallel: parallel kahip')
plt.bar(index + width, sequential_time_parallel_kahip, width, color="#F26F6F", bottom=parallel_time_parallel_kahip, label='sequential: parallel kahip')
plt.bar(index + 2 * width, parallel_time_lpa, width, color="green", label='parallel: lpa')
plt.bar(index + 2 * width, sequential_time_lpa, width, color="#8DEB6E", bottom=parallel_time_lpa, label='sequential: lpa')
plt.xlabel('# Partitions')
plt.ylabel('Run time (s)')
plt.title('Parallel and sequential runtime')
plt.xticks(index + 1.5 * width, partition_sizes)
legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=3)
plt.savefig("graphics/runtime", bbox_extra_artists=(legend,), bbox_inches='tight')
plt.clf()


parallel_size_kahip = [0] * len(partition_sizes)
parallel_size_parallel_kahip = [0] * len(partition_sizes)
parallel_size_lpa = [0] * len(partition_sizes)
for i in range(len(partition_sizes)):
	parallel_size_kahip[i] = graphs[(partition_sizes[i], "kahip")].parallel_size
	parallel_size_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].parallel_size
	parallel_size_lpa[i] = graphs[(partition_sizes[i], "lpa")].parallel_size
index = np.arange(len(partition_sizes))
width = 0.25
plt.bar(index, parallel_size_kahip, width, color="blue", label='kahip')
plt.bar(index + width, parallel_size_parallel_kahip, width, color="red", label='parallel kahip')
plt.bar(index + 2 * width, parallel_size_lpa, width, color="green", label='lpa')
plt.xlabel('# Partitions')
plt.ylabel('Kernel size')
plt.title('Kernel size after parallel reductions')
plt.xticks(index + 1.5 * width, partition_sizes)
plt.legend(loc='upper left')
plt.savefig("graphics/kernel_size")
plt.clf()


