#ifndef PARALLEL_REDUCTIONS_H
#define PARALLEL_REDUCTIONS_H

// #include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "Reduction.h"
#include "SimpleSet.h"
#include "fast_set.h"

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
    parallel_reductions(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> const &vertexPartitions);
    ~parallel_reductions();

    void reduce_graph_parallel();
    void reduce_graph_sequential();

    void ApplyReductions(int const partition, std::vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining, std::vector<int> &tempInt1, std::vector<int> &tempInt2, fast_set &fastSet, double &time, int &isolatedCliqueCount, int &foldedVertexCount, int &removedTwinCount, int &foldedTwinCount, int &removedUnconfinedVerticesCount);
    void UndoReductions(std::vector<Reduction> const &vReductions);
    std::vector<std::vector<int>> getKernel();
    void applyKernelSolution(std::vector<int> kernel_solution);
    void ApplyKernelSolutionToReductions(std::vector<Reduction> const &vReductions);

    std::vector<int> independent_set;

    size_t size() const { return inGraph.Size(); }

    std::vector<SparseArraySet> const& Neighbors()  const { return neighbors;  }

    void SetAllowVertexFolds(bool const allow) { m_bAllowVertexFolds = allow; }

protected: // methods
    bool removeUnconfined(int const partition, int const vertex, ArraySet &remaining, fast_set &closedNeighborhood, std::vector<int> &neighborhood, std::vector<int> &numNeighborsInS, int &removedUnconfinedVerticesCount);
    bool removeTwin(int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, std::vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount);
    bool RemoveIsolatedClique    (int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, std::vector<bool> &vMarkedVertices, int &isolatedCliqueCount);
    bool FoldVertex(int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, int &foldedVertexCount);
    void initReducableVertices(int numPartitions);
    void updateNeighborhood(int const vertex);
    bool isTwoNeighborhoodInSamePartition(int const vertex, int const partition, ArraySet &remaining);
    void updateAllNeighborhoods();

protected: // members
    std::vector<int> graph_to_kernel_map;
    std::vector<int> kernel_solution;
    std::vector<std::vector<std::vector<Reduction>>> AllReductions;
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<SparseArraySet>     neighbors;
    SimpleSet inGraph;
    SimpleSet neighborhoodChanged;
    std::vector<int> partitions;
    std::vector<std::vector<int>> partition_nodes;
    std::vector<ArraySet> inGraphPerPartition;
    SimpleSet boundaryVertices;
#ifdef TIMERS
    clock_t replaceTimer;
    #endif // TIMERS
    bool m_bAllowVertexFolds;
};

#endif //PARALLEL_REDUCTIONS_H
