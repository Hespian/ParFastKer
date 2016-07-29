#include "full_reductions.h"
#include <vector>
#include <memory>
#include "mis_config.h"
#include "parallel_reductions.h"
#include "sequential/branch_and_reduce_algorithm.h"
#include <omp.h>

full_reductions::full_reductions(std::vector<std::vector<int>> &_adj, std::vector<int> _partitions)
: adj(_adj) 
,partitions(_partitions) 
{
	parallel_reducers = std::vector<std::unique_ptr<parallel_reductions>>();
}

void full_reductions::reduce_graph(std::vector<unsigned int> &vertexTimes) {
	std::cout << "Creating object" << std::endl;
	parallel_reducers.push_back(std::unique_ptr<parallel_reductions>(new parallel_reductions(adj, partitions)));
	std::cout << "Finished creating object" << std::endl;
	std::cout << "Before call to parallel reduce_graph" << std::endl;
	parallel_reducers.back()->reduce_graph_parallel(vertexTimes);
	std::cout << "After call to parallel reduce_graph" << std::endl;
	std::cout << "Kernel size after parallel run: " << parallel_reducers.back()->size() << std::endl;
	std::cout << "Before call to sequential reduce_graph" << std::endl;
	parallel_reducers.back()->reduce_graph_sequential();
	std::cout << "After call to sequential reduce_graph" << std::endl;
	std::cout << "Kernel size after sequential run: " << parallel_reducers.back()->size() << std::endl;
}

size_t full_reductions::get_current_is_size_with_folds() {
	// return sequential_reducer->get_current_is_size_with_folds();
	return -1;
}

size_t full_reductions::number_of_nodes_remaining() {
	return parallel_reducers.back()->size();
}

std::vector<std::vector<int>> full_reductions::getKernel() {
	return parallel_reducers.back()->getKernel();
}

void full_reductions::applyKernelSolution(std::vector<int> kernel_solution) {
	parallel_reducers.back()->applyKernelSolution(kernel_solution);
}

std::vector<int> full_reductions::getSolution() {
	return parallel_reducers.front()->independent_set;
}
