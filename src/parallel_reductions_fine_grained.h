 /******************************************************************************
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *****************************************************************************/

#ifndef PARALLEL_REDUCTIONS_FINE_GRAINED_H
#define PARALLEL_REDUCTIONS_FINE_GRAINED_H

// #include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "Reduction.h"
#include "SimpleSet.h"
#include "fast_set.h"
#include "MaximumMatching.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>
#include <string>
#include <atomic>
#include <memory>

#define TIMERS

class parallel_reductions_fine_grained
{
public:
    parallel_reductions_fine_grained(std::vector<std::vector<int>> const &adjacencyArray);
    ~parallel_reductions_fine_grained();

    void reduce_graph_parallel();

    std::vector<std::vector<int>> getKernel();
    void applyKernelSolution(std::vector<int> kernel_solution);
    void ApplyKernelSolutionToReductions(std::vector<Reduction> const &vReductions);

    std::vector<int> independent_set;

    size_t size() const { return inGraph.Size(); }

protected: // methods
    bool RemoveUnconfined(std::vector<fast_set> &closedNeighborhoodPerThread, std::vector<std::vector<int>> &neighborhoodPerThread, std::vector<std::vector<int>> &numNeighborsInSPerThread, std::vector<std::vector<int>> &neighborsInSPerThread, std::vector<char> &isCandidate, std::vector<int> &candidates, std::vector<int> &toRemove, std::vector<std::unique_ptr<std::vector<int>>> &temp);
    bool RemoveIsolatedClique(std::vector<std::vector<bool>> &vMarkedVerticesPerThread, std::vector<int> &toRemove, std::vector<std::unique_ptr<std::vector<int>>> &temp);
    bool LPReduction();

    int degree(int const vertex);

    // Just for testing
    bool checkDegrees();

protected: // members
    std::vector<int> graph_to_kernel_map;
    std::vector<int> kernel_solution;
    std::vector<std::vector<std::vector<Reduction>>> AllReductions;
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<SparseArraySet>     neighbors;
    SimpleSet inGraph;
    MaximumMatching maximumMatching;
    std::vector<std::atomic_int> vertexDegree;
#ifdef TIMERS
    clock_t replaceTimer;
    #endif // TIMERS
    bool m_bAllowVertexFolds;
};

#endif //PARALLEL_REDUCTIONS_FINE_GRAINED_H
