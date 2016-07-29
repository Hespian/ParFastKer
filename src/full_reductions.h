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
	full_reductions(std::vector<std::vector<int>> &_adj, std::vector<int> _partitions);
	void reduce_graph(std::vector<unsigned int> &vertexTimes);
	size_t get_current_is_size_with_folds();
	size_t number_of_nodes_remaining();
	std::vector<std::vector<int>> getKernel();
	void applyKernelSolution(std::vector<int> kernel_solution);
	std::vector<int> getSolution();
private:
	std::vector<std::vector<int>> &adj;
	std::vector<int> partitions;
	std::vector<std::unique_ptr<parallel_reductions>> parallel_reducers;
};


#endif // FULL_REDUCTIONS_H