#ifndef FULL_REDUCTIONS_H
#define FULL_REDUCTIONS_H

#include <vector>
#include <memory>
#include "mis_config.h"
#include "parallel_reductions.h"
#include "sequential/branch_and_reduce_algorithm.h"


class full_reductions
{
public:
	full_reductions(std::vector<std::vector<int>> &_adj, MISConfig &mis_config);
	void reduce_graph();
	size_t get_current_is_size_with_folds();
	size_t number_of_nodes_remaining();
	std::vector<std::vector<int>> getKernel();
	void applyKernelSolution(std::vector<int> kernel_solution);
	std::vector<int> getSolution();
private:
	std::vector<std::vector<int>> &adj;
	MISConfig &mis_config;
	std::vector<std::unique_ptr<parallel_reductions>> parallel_reducers;
	// std::unique_ptr<branch_and_reduce_algorithm> sequential_reducer;
};


#endif // FULL_REDUCTIONS_H