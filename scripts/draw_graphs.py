import matplotlib.pyplot as plt
from collections import namedtuple
import numpy as np

base_path = "../../results/parallel_sequential_58ad3e463/it-2004-sorted.graph"
partition_sizes = [1, 2, 4, 8, 16, 32, 48, 64]
modes = ["parallel_kahip"]

GraphData = namedtuple('GraphData', 'parallel_time parallel_size num_isolated_clique_reductions num_vertex_fold_reductions time_per_partition sequential_time')

def read_file(file):
	file = open(file, "r")
	parallel = True
	for line in file:
		words = line.split()
		if parallel:
			if "|-Nodes:" in line:
				graph_size = int(words[-1])
			if "No. of partitions:" in line:
				num_partitions = int(words[-1])
				num_isolated_clique_reductions = [0] * num_partitions
				num_vertex_fold_reductions = [0] * num_partitions
				time_per_partition = [0.0] * num_partitions
			if "Total time spent applying reductions  : " in line:
				parallel_time = float(words[6])
			if "Kernel size after parallel run: " in line:
				parallel_size = int(words[5])
			if "Currently queued vertices: 0" in line:
				partition = int(words[0][:-1])
				current_partition_num_iloated_cliques = int(words[10][:-1])
				current_partition_num_vertex_folds = int(words[14])
				num_isolated_clique_reductions[partition] = current_partition_num_iloated_cliques
				num_vertex_fold_reductions[partition] = current_partition_num_vertex_folds
			if "Time spent applying reductions" in line:
				partition = int(words[0][:-1])
				time_per_partition[partition] = float(words[-1])
			if "Before call to sequential reduce_graph" in line:
				parallel = False
		else:
			if "Total time spent applying reductions  : " in line:
				sequential_time = float(words[6])

	file.close()

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

	return GraphData(parallel_time, parallel_size, num_isolated_clique_reductions, num_vertex_fold_reductions, time_per_partition, sequential_time)

graphs = {}
for num_partitions in partition_sizes:
	for mode in modes:
		# print(base_path + "-" + str(num_partitions) + "-" + mode)
		graphs[(num_partitions, mode)] = read_file(base_path + "-" + str(num_partitions) + "-" + mode)


sum_vertex_fold_parallel_kahip = [0] * len(partition_sizes)
sum_isolated_clique_parallel_kahip = [0] * len(partition_sizes)
for i in range(len(partition_sizes)):
	sum_vertex_fold_parallel_kahip[i] = sum(graphs[(partition_sizes[i], "parallel_kahip")].num_vertex_fold_reductions)
	sum_isolated_clique_parallel_kahip[i] = sum(graphs[(partition_sizes[i], "parallel_kahip")].num_isolated_clique_reductions)
index = np.arange(len(partition_sizes))
width = 0.75
plt.figure(1)
plt.subplot(211)
plt.bar(index, sum_vertex_fold_parallel_kahip, width, color="red", label='parallel kahip')
plt.xlabel('# Partitions')
plt.ylabel('# Reductions')
plt.title('Number of vertex fold reductions performed in parallel')
plt.subplots_adjust(hspace = 0.5)
plt.subplot(212)
plt.bar(index, sum_isolated_clique_parallel_kahip, width, color="red", label='parallel kahip')
plt.xlabel('# Partitions')
plt.ylabel('# Reductions')
plt.title('Number of isolated clique reductions performed in parallel')
plt.xticks(index, partition_sizes)
plt.savefig("graphics/num_reductions_parallel")
plt.clf()


parallel_time_parallel_kahip = [0] * len(partition_sizes)
sequential_time_parallel_kahip = [0] * len(partition_sizes)
for i in range(len(partition_sizes)):
	parallel_time_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].parallel_time
	sequential_time_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].sequential_time
index = np.arange(len(partition_sizes))
width = 0.75
plt.bar(index, parallel_time_parallel_kahip, width, color="red", label='parallel: parallel kahip')
plt.bar(index, sequential_time_parallel_kahip, width, color="#F26F6F", bottom=parallel_time_parallel_kahip, label='sequential: parallel kahip')
plt.xlabel('# Partitions')
plt.ylabel('Run time (s)')
plt.title('Run time')
plt.xticks(index, partition_sizes)
plt.savefig("graphics/runtime", bbox_inches='tight')
plt.clf()


parallel_size_parallel_kahip = [0] * len(partition_sizes)
for i in range(len(partition_sizes)):
	parallel_size_parallel_kahip[i] = graphs[(partition_sizes[i], "parallel_kahip")].parallel_size
index = np.arange(len(partition_sizes))
width = 0.75
plt.bar(index, parallel_size_parallel_kahip, width, color="red", label='parallel kahip')
plt.xlabel('# Partitions')
plt.ylabel('Kernel size')
plt.title('Kernel size after parallel reductions')
plt.xticks(index, partition_sizes)

plt.savefig("graphics/kernel_size")
plt.clf()


data =  [None] * len(partition_sizes)
for i in range(len(partition_sizes)):
	data[i] = graphs[(partition_sizes[i], "parallel_kahip")].time_per_partition
plt.boxplot(data, labels=partition_sizes, showmeans=True)
plt.xlabel('# Partitions')
plt.ylabel('Run time (s)')
plt.title('Run time per partition box plots\n (Partitioner: Parallel kahip)')
plt.xticks(index, partition_sizes)
plt.savefig("graphics/runtime_boxplots")
plt.clf()


