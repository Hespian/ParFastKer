#include "full_reductions.h"
#include <vector>
#include <memory>
#include "mis_config.h"
#include "parallel_reductions.h"
#include "sequential/branch_and_reduce_algorithm.h"
#include <omp.h>
#include <algorithm>
#include <fstream>

full_reductions::full_reductions(std::vector<std::vector<int>> &_adj, std::vector<int> _partitions)
: adj(_adj) 
,partitions(_partitions) 
{
	parallel_reducers = std::vector<std::unique_ptr<parallel_reductions>>();
}

void full_reductions::reduce_graph() {
    reduce_graph(false, "");
}

void full_reductions::reduce_graph(bool writeKernel, std::string filepath) {
    int numPartitions = *std::max_element(partitions.begin(), partitions.end());
	std::cout << "Creating object" << std::endl;
	parallel_reducers.push_back(std::unique_ptr<parallel_reductions>(new parallel_reductions(adj, partitions)));
	std::cout << "Finished creating object" << std::endl;
	std::cout << "Before call to parallel reduce_graph" << std::endl;
	parallel_reducers.back()->reduce_graph_parallel();
	std::cout << "After call to parallel reduce_graph" << std::endl;
	std::cout << "Kernel size after parallel run: " << parallel_reducers.back()->size() << std::endl;
	std::cout << "Before call to sequential reduce_graph" << std::endl;
    if(writeKernel) {
        writeParallelKernel(filepath);
    }
    if(numPartitions == 0 || numPartitions == 31)
        parallel_reducers.back()->reduce_graph_sequential();
	std::cout << "After call to sequential reduce_graph" << std::endl;
	std::cout << "Kernel size after sequential run: " << parallel_reducers.back()->size() << std::endl;
}

void full_reductions::writeParallelKernel(std::string kernel_path) {
    
    int numKernelEdges = 0;
    std::vector<std::vector<int>> kernel = getKernel();
    for(std::vector<int> neighbors : kernel) {
        numKernelEdges += neighbors.size();
    }
    std::cout << "Writing quasi kernel to " << kernel_path << std::endl;
    std::ofstream f(kernel_path.c_str());
    f << kernel.size() <<  " " <<  numKernelEdges / 2 << std::endl;

    for(auto vertexNeighbors : kernel) {
        std::sort(vertexNeighbors.begin(), vertexNeighbors.end());
        for(auto neighbor : vertexNeighbors) {
            f <<   neighbor + 1 << " " ;
        } 
        f <<  std::endl;
    }

    f.close();

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
