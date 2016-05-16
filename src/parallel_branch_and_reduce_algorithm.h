 /******************************************************************************
 * parallel_branch_and_reduce_algorithm.h
 *
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef PARALLEL_BRANCH_AND_REDUCE_SOLVER_H
#define PARALLEL_BRANCH_AND_REDUCE_SOLVER_H

// local includes
#include "fast_set.h"
#include "parallel_modified.h"
#include "mis_config.h"

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "timer.h"

// system includes
#include <vector>
#include <string>
#include <memory>


class parallel_branch_and_reduce_algorithm
{
friend class parallel_modified;
friend class fold;
friend class alternative;

public:
	static int REDUCTION;
	static int LOWER_BOUND;
	static int BRANCHING;
	static bool outputLP;
	
	std::vector<std::vector<int>> adj;
	static long nBranchings;
	static int debug;
	double SHRINK;
	int    depth;
    int    maxDepth;
    int    rootDepth;
	int    n;
    int    N;


    MISConfig &mis_config;

    /*
     * The partition used for parallelization
     */
     std::vector<int> partitions;
     std::vector<std::vector<int>> partition_nodes;
	/**
	 * current best solution
	 */
    std::vector<int> y;

	/**
	 * current solution (-1: not determined, 0: not in the vc, 1: in the vc, 2: removed by foldings)
	 */
    std::vector<int> x;

	/**
	 * #remaining vertices
	 */
	int rn;
	
	/**
	 * max flow
	 */
	std::vector<int> in;
    std::vector<int> out;
	
	std::vector<std::vector<int>> level;

	std::vector<std::vector<int>> modTmp;

	std::vector<std::vector<std::shared_ptr<parallel_modified>>> modifieds;
	std::vector<int> modifiedN;
			
	parallel_branch_and_reduce_algorithm(std::vector<std::vector<int>> &_adj, int const _N, MISConfig mis_config);

	int deg(int v);
	void set(int v, int a);

	std::vector<fast_set> used;

    // helpers for modifying the graph
	void compute_fold(std::vector<int> const &S, std::vector<int> const &NS, int partition);
	void reverse();

    // reduction methods
	bool fold2Reduction();
	bool isolatedCliqueReduction();
	bool isolatedCliqueReduction(NodeID v, int partition);

    // recursive methods
	bool reduce();

    // vestiges of original Java code
#if 0
	
	void debug(String str, Object...os) {
		StringBuilder sb = new StringBuilder();
		Calendar c = Calendar.getInstance();
		sb.append(String.format("%02d:%02d:%02d  ", c.get(Calendar.HOUR_OF_DAY), c.get(Calendar.MINUTE), c.get(Calendar.SECOND)));
		for (int i = 0; i < depth && i <= maxDepth; i++) sb.append(' ');
		System.err.print(sb);
		System.err.printf(str, os);
	}
#endif // 0
	
    void reduce_graph();

    std::string debugString() const;
    void PrintState() const;

    size_t get_current_is_size() const;
    size_t get_current_is_size_with_folds() const;
    bool   folded_vertices_exist() const;
    std::vector<int> compute_maximal_is();
    size_t compute_alternative_maximal_is_size();
    size_t number_of_nodes_remaining() const;
    void force_into_independent_set(std::vector<NodeID> const &nodes);

#if 0
}
#endif // 0
};

#endif //PARALLEL_BRANCH_AND_REDUCE_SOLVER_H