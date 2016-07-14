#ifndef PARALLEL_REDUCTIONS_H
#define PARALLEL_REDUCTIONS_H

// #include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "Reduction.h"
#include "SimpleSet.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>
#include <string>

#define TIMERS

class parallel_reductions
{
public:
    parallel_reductions(std::vector<std::vector<int>> const &adjacencyArray);
    ~parallel_reductions();

    void reduce_graph(int numPartitions, std::string partitioner);

    void ApplyReductions(int const partition, std::vector<int> vertices, std::vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining);
    void UndoReductions(std::vector<Reduction> const &vReductions);
    std::vector<std::vector<int>> getKernel();
    void applyKernelSolution(std::vector<int> kernel_solution);
    void ApplyKernelSolutionToReductions(std::vector<Reduction> const &vReductions);

    std::vector<int> independent_set;

    size_t size() const { return inGraph.Size(); }

    std::vector<SparseArraySet> const& Neighbors()  const { return neighbors;  }

    void SetAllowVertexFolds(bool const allow) { m_bAllowVertexFolds = allow; }

protected: // methods
    bool RemoveIsolatedClique    (int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, std::vector<bool> &vMarkedVertices, int &isolatedCliqueCount);
    bool FoldVertex(int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, int &foldedVertexCount);
    void partitionGraph(int numPartitions, std::string partitioner);
    void initReducableVertices(int numPartitions);
    void updateNeighborhood(int const vertex);
    bool isTwoNeighborhoodInSamePartition(int const vertex, int const partition);

protected: // members
    std::vector<int> graph_to_kernel_map;
    std::vector<int> kernel_solution;
    std::vector<std::vector<Reduction>> ReductionsPerPartition;
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<SparseArraySet>     neighbors;
    SimpleSet inGraph;
    SimpleSet neighborhoodChanged;
    std::vector<int> partitions;
    std::vector<std::vector<int>> partition_nodes;
#ifdef TIMERS
    clock_t replaceTimer;
    #endif // TIMERS
    bool m_bAllowVertexFolds;
};

#endif //PARALLEL_REDUCTIONS_H
