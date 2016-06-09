#include "full_reductions.h"
#include <vector>
#include <memory>
#include "mis_config.h"
#include "parallel_reductions.h"
#include "sequential/branch_and_reduce_algorithm.h"
#include <omp.h>

full_reductions::full_reductions(std::vector<std::vector<int>> &_adj, MISConfig &mis_config)
: adj(_adj) 
, mis_config(mis_config) {
	parallel_reducers = std::vector<std::unique_ptr<parallel_reductions>>();
}

void full_reductions::reduce_graph() {
	parallel_reducers.push_back(std::unique_ptr<parallel_reductions>(new parallel_reductions(adj, adj.size(), mis_config)));
	parallel_reducers.back()->reduce_graph();	
	std::cout << "Kernel size after parallel run: " << parallel_reducers.back()->number_of_nodes_remaining() << std::endl;
/*	std::vector<std::vector<int> > kernel_adj = parallel_reducers.back()->getKernel();
	parallel_reducers.push_back(std::unique_ptr<parallel_reductions>(new parallel_reductions(kernel_adj, kernel_adj.size(), mis_config)));
	parallel_reducers.back()->reduce_graph();
	std::cout << "Kernel size after second parallel run: " << parallel_reducers.back()->number_of_nodes_remaining() << std::endl;*/
	std::vector<std::vector<int>> kernel_adj = parallel_reducers.back()->getKernel();
	sequential_reducer = std::unique_ptr<branch_and_reduce_algorithm>(new branch_and_reduce_algorithm(kernel_adj, kernel_adj.size()));
	double begin, end, elapsed_secs;
	begin = omp_get_wtime();
	sequential_reducer->reduce_graph_few_reductions();
    end = omp_get_wtime();
    elapsed_secs = double(end - begin);
    std::cout << "Sequential took " << elapsed_secs << " seconds" << std::endl;
}

size_t full_reductions::get_current_is_size_with_folds() {
	return sequential_reducer->get_current_is_size_with_folds();
}

size_t full_reductions::number_of_nodes_remaining() {
	return sequential_reducer->number_of_nodes_remaining();
}

std::vector<std::vector<int>> full_reductions::getKernel() {
	return sequential_reducer->getKernel();
}

void full_reductions::applyKernelSolution(std::vector<int> kernel_solution) {
	sequential_reducer->applyKernelSolution(kernel_solution);
	sequential_reducer->undoReductions();
	std::vector<int> *tmp_kernel_solution = &(sequential_reducer->x);
	for(int i = parallel_reducers.size() - 1; i >= 0; --i) {
		parallel_reducers[i]->applyKernelSolution(*tmp_kernel_solution);
		parallel_reducers[i]->undoReductions();
		tmp_kernel_solution = &(parallel_reducers[i]->x);
	}
}

std::vector<int> full_reductions::getSolution() {
	return parallel_reducers.front()->x;
}
