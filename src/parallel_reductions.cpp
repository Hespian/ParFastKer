#include "parallel_reductions.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "ProfilingHelper.h"

#include <vector>
#include <set>
#include <iostream>
#include <ctime>
#include <cassert>
#include <climits>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "omp.h"
#include <assert.h>

#include "kaHIP_interface.h"

#define ISOLATED_CLIQUE_MAX_NEIGHBORS 3

#define INSERT_REMAINING(partition, remaining, v) if(partitions[v] == partition) remaining.Insert(v);
#define REMOVE_NEIGHBOR(partition, neighbor, vertex) if(partition == partitions[neighbor]) {neighbors[neighbor].Remove(vertex);} else { assert(partition != partitions[neighbor]); neighborhoodChanged.Insert(neighbor);}

using namespace std;

ProfilingHelper_t profilingHelper;

parallel_reductions::parallel_reductions(vector<vector<int>> const &adjacencyArray, vector<int> const &vertexPartitions)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size(), true)
 , neighborhoodChanged(adjacencyArray.size(), false)
 , boundaryVertices(adjacencyArray.size(), false)
 , partitions(vertexPartitions)
#ifdef TIMERS
 , replaceTimer(0)
 #endif // TIMERS
 , m_bAllowVertexFolds(true)
{
    std::cout << "Start constructor" << std::endl;
    int N = adjacencyArray.size();
    for (size_t u=0; u < N; ++u) {
        neighbors[u].InitializeFromAdjacencyArray(m_AdjacencyArray, u);
    }
    independent_set.resize(m_AdjacencyArray.size());
    int numPartitions = *max_element(partitions.begin(), partitions.end()) + 1;
    partition_nodes = std::vector<std::vector<int>>(numPartitions);
    for(int node = 0; node < N; ++node) {
        assert(partitions[node] >= 0);
        assert(partitions[node] < numPartitions);
        partition_nodes[partitions[node]].push_back(node);
    }

    for(int node = 0; node < N; ++node) {
        boundaryVertices.Remove(node);
        for(auto neighbor: neighbors[node]) {
            if(partitions[neighbor] != partitions[node]) {
                boundaryVertices.Insert(node);
            }
        }
    }
    for(int partition = 0; partition < numPartitions; ++partition) {
        std::cout << partition << ": " << partition_nodes[partition].size() << " vertices" << std::endl;
    }
    std::cout << "Finished constructor" << std::endl;
}

parallel_reductions::~parallel_reductions()
{

#ifdef TIMERS
    cout << "Total time spent undoing  reductions  : " << (replaceTimer/(double)CLOCKS_PER_SEC) << endl;
#endif // TIMERS
}

std::vector<std::vector<int>> parallel_reductions::getKernel() {
    graph_to_kernel_map = std::vector<int> (m_AdjacencyArray.size());
    int nodecount = 0;
    for(int node = 0; node < m_AdjacencyArray.size(); node++) {
        if(inGraph.Contains(node)) {
            graph_to_kernel_map[node] = nodecount++;
        }
    }

    std::vector<std::vector<int>> kernel_adj(nodecount);

    // Build adjacency vectors
    for(int node = 0; node < m_AdjacencyArray.size(); node++) {
        if(inGraph.Contains(node)) {
            kernel_adj[graph_to_kernel_map[node]].reserve(neighbors[node].Size());
            for(auto neighbor : neighbors[node]) {
                if(inGraph.Contains(neighbor)) {
                    kernel_adj[graph_to_kernel_map[node]].push_back(graph_to_kernel_map[neighbor]);
                }
            }
            std::sort(kernel_adj[graph_to_kernel_map[node]].begin(), kernel_adj[graph_to_kernel_map[node]].end());
        }
    }
    return kernel_adj;
}

void parallel_reductions::applyKernelSolution(std::vector<int> kernel_solution){
    for(int node = 0; node < m_AdjacencyArray.size(); ++node) {
        if(inGraph.Contains(node)) {
            independent_set[node] = kernel_solution[graph_to_kernel_map[node]];
        }
    }
    for(int i = AllReductions.size(); i > 0; i--) {
        for(auto Reductions: AllReductions[i - 1]) {
            ApplyKernelSolutionToReductions(Reductions);
        }
    }
}

bool parallel_reductions::RemoveIsolatedClique(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &isolatedCliqueCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    if(boundaryVertices.Contains(vertex)) {
        profilingAddTimeUnsuccessfulIsolatedCliquePartition(&profilingHelper, partition);
        return false;
    }

    /*if(neighbors[vertex].Size() > ISOLATED_CLIQUE_MAX_NEIGHBORS)
        return false;*/

    for (int const neighbor : neighbors[vertex]) {
        assert(partitions[neighbor] == partition);
        if (neighbors[neighbor].Size() < neighbors[vertex].Size()) {
            profilingAddTimeUnsuccessfulIsolatedCliqueDegree(&profilingHelper, partition);
            return false;
        }
    }

    bool superSet(true);

    for (int const neighbor : neighbors[vertex]) {

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = true;
        }
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) {
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        for (int const nNeighbor : neighbors[neighbor]) {
            vMarkedVertices[nNeighbor] = false;
        }
        vMarkedVertices[neighbor] = false;

        if (!superSet) {
            profilingAddTimeUnsuccessfulIsolatedCliqueNoClique(&profilingHelper, partition);
            return false;
        }
    }
    if (superSet) {
        // Reduction reduction(ISOLATED_VERTEX);
        // reduction.SetVertex(vertex);
        independent_set[vertex] = 0;
        for (int const neighbor : neighbors[vertex]) {
            inGraph.Remove(neighbor);
            boundaryVertices.Remove(neighbor);
            remaining.Remove(neighbor);
            independent_set[neighbor] = 1;
            // reduction.AddNeighbor(neighbor);
            // reduction.AddRemovedEdge(vertex,   neighbor);
            // reduction.AddRemovedEdge(neighbor, vertex);
        }
        inGraph.Remove(vertex);
        boundaryVertices.Remove(vertex);

        for (int const neighbor : neighbors[vertex]) {
            for (int const nNeighbor : neighbors[neighbor]) {
                if (inGraph.Contains(nNeighbor)) {
                    INSERT_REMAINING(partition, remaining, nNeighbor);
                    // remaining.Insert(nNeighbor);
                }

                if (nNeighbor != vertex) {
                    REMOVE_NEIGHBOR(partition, nNeighbor, neighbor);
                    // neighbors[nNeighbor].Remove(neighbor);
                    // reduction.AddRemovedEdge(nNeighbor, neighbor);
                    // reduction.AddRemovedEdge(neighbor, nNeighbor);
                }
            }
            neighbors[neighbor].Clear();
        }
        neighbors[vertex].Clear();

        // vReductions.emplace_back(std::move(reduction));
        isolatedCliqueCount++;

        profilingAddTimeSuccessfulIsolatedClique(&profilingHelper, partition);
        return true;
    }
    assert(false);
    return false;
}

bool parallel_reductions::isTwoNeighborhoodInSamePartition(int const vertex, int const partition, ArraySet &remaining) {
    if(boundaryVertices.Contains(vertex)) {
        return false;
    }
    for(int neighbor : neighbors[vertex]) {
        bool wasBoundaryVertex = boundaryVertices.Contains(neighbor);
        updateNeighborhood(neighbor);
        if(boundaryVertices.Contains(neighbor)) {
            return false;
        } else if(wasBoundaryVertex) {
            remaining.Insert(neighbor);
        }
    }
    return true;
}

bool parallel_reductions::FoldVertex(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, int &foldedVertexCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    if (neighbors[vertex].Size() != 2) { 
        profilingAddTimeUnsuccessfulFoldDegree(&profilingHelper, partition);
        return false;
    }
    if (!isTwoNeighborhoodInSamePartition(vertex, partition, remaining)) { 
        profilingAddTimeUnsuccessfulFoldWrongPartition(&profilingHelper, partition);
        return false;
    }
    if (neighbors[neighbors[vertex][0]].Contains(neighbors[vertex][1])) {
        profilingAddTimeUnsuccessfulFoldAdjacent(&profilingHelper, partition);
        return false; // neighbors can't be adjacent.
    }

    foldedVertexCount++;

    int const vertex1(neighbors[vertex][0]);
    int const vertex2(neighbors[vertex][1]);

    assert(partitions[vertex1] == partition);
    assert(partitions[vertex2] == partition);

    Reduction reduction(FOLDED_VERTEX);
    reduction.SetVertex(vertex);
    reduction.AddNeighbor(vertex1);
    reduction.AddNeighbor(vertex2);

    neighbors[vertex].Clear();
    neighbors[vertex].Resize(neighbors[vertex1].Size() + neighbors[vertex2].Size());
    neighbors[vertex1].Remove(vertex);
    neighbors[vertex2].Remove(vertex);

    // reduction.AddRemovedEdge(vertex, vertex1);
    // reduction.AddRemovedEdge(vertex1, vertex);
    // reduction.AddRemovedEdge(vertex, vertex2);
    // reduction.AddRemovedEdge(vertex2, vertex);

    for (int const neighbor1 : neighbors[vertex1]) {
        assert(partitions[neighbor1] == partition);
        if (neighbor1 == vertex) continue;
        neighbors[neighbor1].Remove(vertex1);
        // reduction.AddRemovedEdge(neighbor1, vertex1);
        // reduction.AddRemovedEdge(vertex1, neighbor1);
        neighbors[vertex].Insert(neighbor1);
        assert(partitions[vertex] == partitions[neighbor1]);
        INSERT_REMAINING(partition, remaining, neighbor1);
        // remaining.Insert(neighbor1);
    }
    neighbors[vertex1].Clear();
    assert(partitions[vertex] == partitions[vertex1]);

    for (int const neighbor2 : neighbors[vertex2]) {
        assert(partitions[neighbor2] == partition);
        if (neighbor2 == vertex) continue;
        neighbors[neighbor2].Remove(vertex2);
        // reduction.AddRemovedEdge(neighbor2, vertex2);
        // reduction.AddRemovedEdge(vertex2, neighbor2);
        neighbors[vertex].Insert(neighbor2);
        assert(partitions[vertex] == partitions[neighbor2]);
        INSERT_REMAINING(partition, remaining, neighbor2);
        // remaining.Insert(neighbor2);
    }

    neighbors[vertex2].Clear();
    assert(partitions[vertex] == partitions[vertex2]);
    
    for (int const neighbor : neighbors[vertex]) {
        neighbors[neighbor].Insert(vertex);
    }

    INSERT_REMAINING(partition, remaining, vertex);
    // remaining.Insert(vertex);

    vReductions.emplace_back(std::move(reduction));

    remaining.Remove(vertex1);
    remaining.Remove(vertex2);
    inGraph.Remove(vertex2);
    inGraph.Remove(vertex1);

    profilingAddTimeSuccessfulFold(&profilingHelper, partition);
    return true;
}

void parallel_reductions::updateNeighborhood(int const vertex) {
    if(!neighborhoodChanged.Contains(vertex))
        return;
    neighborhoodChanged.Remove(vertex);

    profilingStartClockUpdateNeighborhood(&profilingHelper, partitions[vertex], vertex);

    std::vector<int> verticesToRemove;
    int partition = partitions[vertex];
    bool isBoundaryVertex = false;
    for(int neighbor: neighbors[vertex]) {
        if(!inGraph.Contains(neighbor)) {
            verticesToRemove.push_back(neighbor);
        } else if(partitions[neighbor] != partition) {
            isBoundaryVertex = true;
        }
    }
    for(int neighborToRemove : verticesToRemove) {
        neighbors[vertex].Remove(neighborToRemove);
    }
    if(!isBoundaryVertex) {
        boundaryVertices.Remove(vertex);
    }
    profilingAddTimeUpdateNeighborhood(&profilingHelper, partitions[vertex]);
}

void parallel_reductions::reduce_graph_parallel() {
    int numPartitions = partition_nodes.size();
    profilingInit(&profilingHelper, &neighbors, numPartitions);

    vector<vector<bool>> vMarkedVerticesPerPartition(numPartitions);
    vector<ArraySet> remainingPerPartition(numPartitions);
    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(numPartitions);
    for(int partition = 0; partition < numPartitions; partition++) {
        std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
        vMarkedVerticesPerPartition[partition] = vMarkedVertices;
        ArraySet remaining = ArraySet(m_AdjacencyArray.size());
        remainingPerPartition[partition] = remaining;
    }

    vector<double> partitionTimes(numPartitions);
    
    double startClock = omp_get_wtime();

    #pragma omp parallel for
    for(int partition = 0; partition < numPartitions; partition++) {
        std::cout << partition << ": Filling remaining vertices set..." << std::endl;
        remainingPerPartition[partition].Clear();
        for (int const vertex : partition_nodes[partition]) {
            if(inGraph.Contains(vertex)) {
                assert(partitions[vertex] == partition);
                remainingPerPartition[partition].Insert(vertex);
            }
        }
        ApplyReductions(partition, ReductionsPerPartition[partition], vMarkedVerticesPerPartition[partition], remainingPerPartition[partition], partitionTimes[partition]);
    }


    double endClock = omp_get_wtime();
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Time spent applying reductions  : " << partitionTimes[partition] << endl;
    }
    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
}

void parallel_reductions::reduce_graph_sequential() {
    profilingInit(&profilingHelper, &neighbors, 1);

    int N = m_AdjacencyArray.size();

    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(1);
    std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
    ArraySet remaining = ArraySet(N);

    partition_nodes.resize(1);
    partition_nodes[0].clear();
    for(int vertex = 0; vertex < N; ++vertex) {
        boundaryVertices.Remove(vertex);
        partitions[vertex] = 0;
        partition_nodes[0].push_back(vertex);
    }

    updateAllNeighborhoods();

    double time(0);
    
    double startClock = omp_get_wtime();

    std::cout << "Filling remaining vertices set..." << std::endl;
    remaining.Clear();
    for(int vertex = 0; vertex < N; ++vertex)  {
        if(inGraph.Contains(vertex)) {
            assert(partitions[vertex] == 0);
            remaining.Insert(vertex);
        }
    }
    ApplyReductions(0, ReductionsPerPartition[0], vMarkedVertices, remaining, time);


    double endClock = omp_get_wtime();
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
}

void parallel_reductions::updateAllNeighborhoods() {
    double startClock = omp_get_wtime();
    #pragma omp parallel for
    for(int vertex = 0; vertex < neighbors.size(); ++vertex)
        updateNeighborhood(vertex);
    double endClock = omp_get_wtime();
    cout << "Time spent updating neighborhoods  : " << (endClock - startClock) << endl;
}

void parallel_reductions::ApplyReductions(int const partition, vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining, double &time)
{
    double startClock = omp_get_wtime();
    int iterations(0);
    int foldedVertexCount(0);
    int isolatedCliqueCount(0);
    std::cout << partition << ": Starting reductions..." << std::endl;
    while (!remaining.Empty()) {
        int const vertex = *(remaining.begin());
        remaining.Remove(vertex);
        assert(inGraph.Contains(vertex));
        assert(partitions[vertex] == partition);
        updateNeighborhood(vertex);
        bool reduction = RemoveIsolatedClique(partition, vertex, vReductions, remaining, vMarkedVertices, isolatedCliqueCount);
        if (!reduction && m_bAllowVertexFolds) {
            reduction = FoldVertex(partition, vertex, vReductions, remaining, foldedVertexCount);
        }
        iterations++;
        if(iterations % 1000000 == 0) {
            std::cout << partition << ": " << iterations << " iterations. Currently queued vertices: " << remaining.Size() << ". Isolated clique reductions: " << isolatedCliqueCount << ", vertex fold count: " << foldedVertexCount << std::endl;
        }
    }
    std::cout << partition << ": " << iterations << " iterations. Currently queued vertices: " << remaining.Size() << ". Isolated clique reductions: " << isolatedCliqueCount << ", vertex fold count: " << foldedVertexCount << std::endl;
    std::cout << partition << ": Finished reductions!" << std::endl;
    double endClock = omp_get_wtime();
    time = (endClock - startClock);
}

void parallel_reductions::UndoReductions(vector<Reduction> const &vReductions)
{
#ifdef TIMERS
    clock_t startClock = clock();
#endif //TIMERS
    for (size_t index = vReductions.size(); index > 0; index--) {
        Reduction const &reduction(vReductions[index-1]);
        switch(reduction.GetType()) {
            case ISOLATED_VERTEX:

                inGraph.Insert(reduction.GetVertex());
                independent_set[reduction.GetVertex()] = 0;
                for (int const neighbor : reduction.GetNeighbors()) {
                    inGraph.Insert(neighbor);
                    independent_set[neighbor] = 1;
                }
                for (pair<int,int> const &edge : reduction.GetRemovedEdges()) {
                    neighbors[edge.first].Insert(edge.second);
                }
            break;
            case FOLDED_VERTEX:
                inGraph.Insert(reduction.GetNeighbors()[0]);
                inGraph.Insert(reduction.GetNeighbors()[1]);
                // first remove all added edges
                for (int const neighbor : neighbors[reduction.GetVertex()]) {
                    neighbors[neighbor].Remove(reduction.GetVertex());
                }
                neighbors[reduction.GetVertex()].Clear();
                // then replace all removed edges
                for (pair<int,int> const &edge : reduction.GetRemovedEdges()) {
                    neighbors[edge.first].Insert(edge.second);
                }
                if(independent_set[reduction.GetVertex()] == 0) {
                    independent_set[reduction.GetNeighbors()[0]] = 0;
                    independent_set[reduction.GetNeighbors()[1]] = 0;
                    independent_set[reduction.GetVertex()] = 1;
                } else {
                    independent_set[reduction.GetNeighbors()[0]] = 1;
                    independent_set[reduction.GetNeighbors()[1]] = 1;
                    independent_set[reduction.GetVertex()] = 0;
                }
            break;
            default:
                cout << "ERROR!: Unsupported reduction type used..." << endl << flush;
                exit(1);
            break;
        };
    }

#ifdef TIMERS
    clock_t endClock = clock();
    replaceTimer += (endClock - startClock);
    #endif // TIMERS
}

void parallel_reductions::ApplyKernelSolutionToReductions(vector<Reduction> const &vReductions)
{
#ifdef TIMERS
    clock_t startClock = clock();
#endif //TIMERS
    for (size_t index = vReductions.size(); index > 0; index--) {
        Reduction const &reduction(vReductions[index-1]);
        switch(reduction.GetType()) {
            /*case ISOLATED_VERTEX:
                independent_set[reduction.GetVertex()] = 0;
                for (int const neighbor : reduction.GetNeighbors()) {
                    independent_set[neighbor] = 1;
                }
            break;*/
            case FOLDED_VERTEX:
                if(independent_set[reduction.GetVertex()] == 0) {
                    independent_set[reduction.GetNeighbors()[0]] = 0;
                    independent_set[reduction.GetNeighbors()[1]] = 0;
                    independent_set[reduction.GetVertex()] = 1;
                } else {
                    independent_set[reduction.GetNeighbors()[0]] = 1;
                    independent_set[reduction.GetNeighbors()[1]] = 1;
                    independent_set[reduction.GetVertex()] = 0;
                }
            break;
            default:
                cout << "ERROR!: Unsupported reduction type used..." << endl << flush;
                exit(1);
            break;
        };
    }

#ifdef TIMERS
    clock_t endClock = clock();
    replaceTimer += (endClock - startClock);
    #endif // TIMERS
}
