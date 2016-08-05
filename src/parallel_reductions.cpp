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
#include <limits.h>
#include <functional>

#include "kaHIP_interface.h"

#define ISOLATED_CLIQUE_MAX_NEIGHBORS 3

#define INSERT_REMAINING(partition, remaining, v) if(partitions[v] == partition) remaining.Insert(v);
// Remove vertex from inGraph first!
#define REMOVE_NEIGHBOR(partition, neighbor, vertex) {if(partition == partitions[neighbor]) {neighbors[neighbor].Remove(vertex);} else { assert(partition != partitions[neighbor]); neighborhoodChanged.Insert(neighbor);}}

using namespace std;

ProfilingHelper_t profilingHelper;

parallel_reductions::parallel_reductions(vector<vector<int>> const &adjacencyArray, vector<int> const &vertexPartitions)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size(), true)
 , neighborhoodChanged(adjacencyArray.size(), false)
 , boundaryVertices(adjacencyArray.size(), false)
 , partitions(vertexPartitions)
 , independent_set(adjacencyArray.size(), -1)
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
            assert(independent_set[node] == -1);
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

bool parallel_reductions::removeTwin(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount)
{
    assert(partitions[vertex] == partition);
    assert(vMarkedVertices.size() == neighbors.size());
    // This takes really long (it's O(n))
    // assert(std::accumulate(vMarkedVertices.begin(), vMarkedVertices.end(), false, std::logical_or<bool>()) == false);
    if(boundaryVertices.Contains(vertex))
        return false;

    if(neighbors[vertex].Size() != 3)
        return false;

    vector<int> markedVertices;

    int smallestDegreeNeighbor = -1;
    int smallesDegreeNeighborDegree = INT_MAX;
    for(int neighbor: neighbors[vertex]) {
        assert(partitions[neighbor] == partition);
        assert(neighbor != vertex);
        updateNeighborhood(neighbor);
        vMarkedVertices[neighbor] = true;
        markedVertices.push_back(neighbor);
        int neighborDegree = neighbors[neighbor].Size();
        if(neighborDegree < smallesDegreeNeighborDegree) {
            smallesDegreeNeighborDegree = neighborDegree;
            smallestDegreeNeighbor = neighbor;
        }
    }
    assert(smallestDegreeNeighbor != -1);

    int twin = -1;
    for(int possibleTwin: neighbors[smallestDegreeNeighbor]) {
        if(possibleTwin == vertex) continue;
        if(partitions[possibleTwin] != partition) continue;
        updateNeighborhood(possibleTwin);
        if(neighbors[possibleTwin].Size() != 3) continue;
        if(vMarkedVertices[possibleTwin]) continue;
        assert(partitions[possibleTwin] == partitions[vertex]);
        bool isTwin = true;
        for(int twinNeighbor: neighbors[possibleTwin]) {
            if(!vMarkedVertices[twinNeighbor]) {
                isTwin = false;
                break;
            }
        }
        if(isTwin) {
            twin = possibleTwin;
            break;
        }
    }
    if(twin == -1) {
        for(int markedVertex: markedVertices) {
            vMarkedVertices[markedVertex] = false;
        }
        return false;
    }
    assert(twin >= 0);
    assert(neighbors[vertex].Equals(neighbors[twin]));
    assert(partitions[twin] == partitions[vertex]);
    assert(!boundaryVertices.Contains(twin));
    assert(!neighbors[vertex].Contains(twin));

    bool isNeighborhoodAdjacent = false;
    for(int neighbor1: neighbors[vertex]) {
        for(int neighbor2: neighbors[neighbor1]) {
            if(vMarkedVertices[neighbor2]) {
                isNeighborhoodAdjacent = true;
                goto afterNeighborhoodCheck;
            }
        }
    }
afterNeighborhoodCheck:

    bool reduced = false;
    if(isNeighborhoodAdjacent) {
        // Case where all vertices get removed from the graph and the twins get inserted into the independent set
        for(int neighbor1: neighbors[vertex]) {
            assert(partitions[neighbor1] == partition);
            inGraph.Remove(neighbor1);
            for(int neighbor2: neighbors[neighbor1]) {
                if(neighbor2 != vertex && inGraph.Contains(neighbor2)) {
                    REMOVE_NEIGHBOR(partition, neighbor2, neighbor1);
                    INSERT_REMAINING(partition, remaining, neighbor2);
                }
            }
            neighbors[neighbor1].Clear();
            boundaryVertices.Remove(neighbor1);
            remaining.Remove(neighbor1);
            independent_set[neighbor1] = 1;
        }
        neighbors[vertex].Clear();
        inGraph.Remove(vertex);
        independent_set[vertex] = 0;
        remaining.Remove(vertex);
        neighbors[twin].Clear();
        inGraph.Remove(twin);
        independent_set[twin] = 0;
        remaining.Remove(twin);
        removedTwinCount++;
        reduced = true;
    } else {
        // Case where the vertices get folded
        bool twoNeighborHoodInSamePartition = isTwoNeighborhoodInSamePartition(vertex, partition, remaining);
        if(!twoNeighborHoodInSamePartition) {
            reduced = false;
        } else {
            Reduction reduction(FOLDED_TWINS);
            reduction.SetVertex(vertex);
            reduction.SetTwin(twin);
            for(int neighbor: neighbors[vertex])
                reduction.AddNeighbor(neighbor);

            int neighborHoodSize(0);
            for(int neighbor: reduction.GetNeighbors()) {
                assert(partitions[neighbor] == partitions[twin]);
                neighbors[neighbor].Remove(twin);
                neighbors[neighbor].Remove(vertex);
                neighborHoodSize += neighbors[neighbor].Size();
            }
            neighbors[twin].Clear();
            neighbors[vertex].Clear();
            neighbors[vertex].Resize(neighborHoodSize);
            for(int neighbor1: reduction.GetNeighbors()) {
                assert(!boundaryVertices.Contains(neighbor1));
                for(int neighbor2: neighbors[neighbor1]) {
                    assert(partitions[neighbor2] == partitions[neighbor1]);
                    assert(neighbor2 != vertex);
                    assert(neighbor2 != twin);
                    assert(!vMarkedVertices[neighbor2]);
                    neighbors[neighbor2].Remove(neighbor1);
                    neighbors[vertex].Insert(neighbor2);
                    neighbors[neighbor2].Insert(vertex);
                    remaining.Insert(neighbor2);
                }
                neighbors[neighbor1].Clear();
                inGraph.Remove(neighbor1);
                boundaryVertices.Remove(neighbor1);
                remaining.Remove(neighbor1);
            }
            inGraph.Remove(twin);
            assert(!boundaryVertices.Contains(twin));
            remaining.Remove(twin);
            remaining.Insert(vertex);
            vReductions.push_back(reduction);
            assert(inGraph.Contains(vertex));
            assert(!inGraph.Contains(reduction.GetNeighbors()[0]));
            assert(!inGraph.Contains(reduction.GetNeighbors()[1]));
            assert(!inGraph.Contains(reduction.GetNeighbors()[2]));
            foldedTwinCount++;
            reduced = true;
        }
    }

    for(int markedVertex: markedVertices) {
        vMarkedVertices[markedVertex] = false;
    }
    return reduced;
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
    if(!neighborhoodChanged.Contains(vertex)) {
        return;
    }
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
    vector<int> numIsolatedCliqueReductions(numPartitions, 0);
    vector<int> numVertexFoldReductions(numPartitions, 0);
    vector<int> numTwinReductionsRemoved(numPartitions, 0);
    vector<int> numTwinReductionsFolded(numPartitions, 0);
    
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
        ApplyReductions(partition, ReductionsPerPartition[partition], vMarkedVerticesPerPartition[partition], remainingPerPartition[partition], partitionTimes[partition], numIsolatedCliqueReductions[partition], numVertexFoldReductions[partition], numTwinReductionsRemoved[partition], numTwinReductionsFolded[partition]);
    }


    double endClock = omp_get_wtime();
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Time spent applying reductions  : " << partitionTimes[partition] << endl;
    }

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Number of isolated clique reductions: " << numIsolatedCliqueReductions[partition]<< endl;
    }

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Number of vertex fold reductions: " << numVertexFoldReductions[partition]<< endl;
    }

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Number of twin reductions (removed): " << numTwinReductionsRemoved[partition]<< endl;
    }

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Number of twin reductions (folded): " << numTwinReductionsFolded[partition]<< endl;
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
    int numIsolatedCliqueReductions(0);
    int numVertexFoldReductions(0);
    int numTwinReductionsRemoved(0);
    int numTwinReductionsFolded(0);
    
    double startClock = omp_get_wtime();

    std::cout << "Filling remaining vertices set..." << std::endl;
    remaining.Clear();
    for(int vertex = 0; vertex < N; ++vertex)  {
        if(inGraph.Contains(vertex)) {
            assert(partitions[vertex] == 0);
            remaining.Insert(vertex);
        }
    }
    ApplyReductions(0, ReductionsPerPartition[0], vMarkedVertices, remaining, time, numIsolatedCliqueReductions, numVertexFoldReductions, numTwinReductionsRemoved, numTwinReductionsFolded);


    double endClock = omp_get_wtime();
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
    cout << "Number of isolated clique reductions: " << numIsolatedCliqueReductions << endl;
    cout << "Number of vertex fold reductions: " << numVertexFoldReductions << endl;
    cout << "Number of twin reductions (removed): " << numTwinReductionsRemoved << endl;
    cout << "Number of twin reductions (folded): " << numTwinReductionsFolded << endl;
}

void parallel_reductions::updateAllNeighborhoods() {
    double startClock = omp_get_wtime();
    #pragma omp parallel for
    for(int vertex = 0; vertex < neighbors.size(); ++vertex)
        updateNeighborhood(vertex);
    double endClock = omp_get_wtime();
    cout << "Time spent updating neighborhoods  : " << (endClock - startClock) << endl;
}

void parallel_reductions::ApplyReductions(int const partition, vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining, double &time, int &isolatedCliqueCount, int &foldedVertexCount, int &removedTwinCount, int &foldedTwinCount)
{
    double startClock = omp_get_wtime();
    int iterations(0);
    std::cout << partition << ": Starting reductions..." << std::endl;
    while (!remaining.Empty()) {
        int const vertex = *(remaining.begin());
        remaining.Remove(vertex);
        assert(inGraph.Contains(vertex));
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        updateNeighborhood(vertex);
        bool reduction = RemoveIsolatedClique(partition, vertex, vReductions, remaining, vMarkedVertices, isolatedCliqueCount);
        if (!reduction && m_bAllowVertexFolds) {
            reduction = FoldVertex(partition, vertex, vReductions, remaining, foldedVertexCount);
        }
        if (!reduction) {
            reduction = removeTwin(partition, vertex, vReductions, remaining, vMarkedVertices, removedTwinCount, foldedTwinCount);
        }
        iterations++;
        if(iterations % 1000000 == 0) {
            std::cout << partition << ": " << iterations << " iterations. Currently queued vertices: " << remaining.Size() << ". Isolated clique reductions: " << isolatedCliqueCount << ", vertex fold count: " << foldedVertexCount << ", twin reduction count (removed): " << removedTwinCount  << ", twin reduction coint (folded): " << foldedTwinCount << std::endl;
        }
    }
    std::cout << partition << ": " << iterations << " iterations. Currently queued vertices: " << remaining.Size() << ". Isolated clique reductions: " << isolatedCliqueCount << ", vertex fold count: " << foldedVertexCount << ", twin reduction count: " << removedTwinCount  << ", twin reduction coint (folded): " << foldedTwinCount  << std::endl;
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
                assert(independent_set[reduction.GetNeighbors()[0]] == -1);
                assert(independent_set[reduction.GetNeighbors()[1]] == -1);
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
            case FOLDED_TWINS:
                assert(reduction.GetNeighbors().size() == 3);
                assert(independent_set[reduction.GetNeighbors()[0]] == -1);
                assert(independent_set[reduction.GetNeighbors()[1]] == -1);
                assert(independent_set[reduction.GetNeighbors()[2]] == -1);
                assert(independent_set[reduction.GetTwin()] == -1);
                if(independent_set[reduction.GetVertex()] == 0) {
                    independent_set[reduction.GetVertex()] = 1;
                    independent_set[reduction.GetTwin()] = 1;
                    independent_set[reduction.GetNeighbors()[0]] = 0;
                    independent_set[reduction.GetNeighbors()[1]] = 0;
                    independent_set[reduction.GetNeighbors()[2]] = 0;
                } else {
                    independent_set[reduction.GetVertex()] = 0;
                    independent_set[reduction.GetTwin()] = 0;
                    independent_set[reduction.GetNeighbors()[0]] = 1;
                    independent_set[reduction.GetNeighbors()[1]] = 1;
                    independent_set[reduction.GetNeighbors()[2]] = 1;
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
