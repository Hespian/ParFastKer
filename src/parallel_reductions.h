#ifndef PARALLEL_REDUCTIONS_H
#define PARALLEL_REDUCTIONS_H

// #include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "Reduction.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>

#define TIMERS

class parallel_reductions
{
public:
    parallel_reductions(std::vector<std::vector<int>> const &adjacencyArray);
    ~parallel_reductions();

    void reduce_graph();

    void ApplyReductions(std::vector<Reduction> &vReductions);
    void UndoReductions(std::vector<Reduction> const &vReductions);
    std::vector<std::vector<int>> getKernel();
    void applyKernelSolution(std::vector<int> kernel_solution);
    void ApplyKernelSolutionToReductions(std::vector<Reduction> const &vReductions);

    std::vector<int> x;

    size_t size() const { return inGraph.Size(); }

    ArraySet const& GetInGraph()  const { return inGraph;  }
    std::vector<SparseArraySet> const& Neighbors()  const { return neighbors;  }

    size_t GetFoldedVertexCount() const { return foldedVertexCount; }

    void SetAllowVertexFolds(bool const allow) { m_bAllowVertexFolds = allow; }

protected: // methods
    bool RemoveIsolatedClique    (int const vertex, std::vector<Reduction> &vReductions);
    bool FoldVertex(int const vertex, std::vector<Reduction> &vReductions);

protected: // members
    std::vector<int> graph_to_kernel_map;
    std::vector<int> kernel_solution;
    std::vector<Reduction> Reductions;
    std::vector<std::vector<int>> const &m_AdjacencyArray;
    std::vector<SparseArraySet>     neighbors;
    ArraySet inGraph;
    ArraySet remaining;
    std::vector<bool> vMarkedVertices;
#ifdef TIMERS
    clock_t removeTimer;
    clock_t replaceTimer;
    #endif // TIMERS
    size_t foldedVertexCount;
    bool m_bAllowVertexFolds;
};

#endif //PARALLEL_REDUCTIONS_H
