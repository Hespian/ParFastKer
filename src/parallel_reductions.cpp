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

#define ISOLATED_CLIQUE_MAX_NEIGHBORS 2
#define MAX_SIZE_UNCONFINED 6
#define DEPENDENCY_CHECKING_THRESHOLD_MULTIPLIER 3.0
#define DEPENDENCYCHECKING_BURST_ESTIMATION_ALPHA 0.5
#define GLOBAL_BURST_ESTIMATION_ALPHA 0.5
#define GLOBAL_BURST_THRESHOLD_MULTIPLIER 6

#define INSERT_REMAINING(partition, remaining, v) if(partitions[v] == partition) remaining.Insert(v);
// Remove vertex from inGraph first!
#define REMOVE_NEIGHBOR(partition, neighbor, vertex) {assert(!inGraph.Contains(vertex)); assert(partitions[vertex] == partition); vertexDegree[neighbor]--; if(partition != partitions[neighbor]) {numCutEdges[neighbor]--; neighborhoodChanged.Insert(neighbor);}}
#define REMOVE_VERTEX(partition, vertex) {inGraph.Remove(vertex); inGraphPerPartition[partition].Remove(vertex);}

using namespace std;

ProfilingHelper_t profilingHelper;

parallel_reductions::parallel_reductions(vector<vector<int>> const &adjacencyArray, vector<int> const &vertexPartitions)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size(), true)
 , neighborhoodChanged(adjacencyArray.size(), false)
 , partitions(vertexPartitions)
 , independent_set(adjacencyArray.size(), -1)
 , maximumMatching(adjacencyArray)
 , vertexDegree(adjacencyArray.size())
 , numCutEdges(adjacencyArray.size())
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
        numCutEdges[node] = 0;
    }

    for(int node = 0; node < N; ++node) {
        vertexDegree[node] = neighbors[node].Size();
        for(auto neighbor: neighbors[node]) {
            if(partitions[neighbor] != partitions[node]) {
                numCutEdges[node]++;
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

int parallel_reductions::degree(int const vertex) {
    return vertexDegree[vertex];
}

bool parallel_reductions::isBoundaryVertex(const int vertex) {
    return numCutEdges[vertex] > 0;
}

bool parallel_reductions::RemoveIsolatedClique(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &isolatedCliqueCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    if(isBoundaryVertex(vertex)) {
        profilingAddTimeUnsuccessfulIsolatedCliquePartition(&profilingHelper, partition);
        return false;
    }

    int const deg = degree(vertex);

    if(deg > ISOLATED_CLIQUE_MAX_NEIGHBORS)
        return false;

    for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        assert(partitions[neighbor] == partition);
        if (degree(neighbor) < deg) {
            profilingAddTimeUnsuccessfulIsolatedCliqueDegree(&profilingHelper, partition);
            return false;
        }
    }

    bool superSet(true);

    for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        // TODO Should be possible faster
        for (int const nNeighbor : neighbors[neighbor]) if(inGraph.Contains(nNeighbor)) {
            vMarkedVertices[nNeighbor] = true;
        }
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) if(inGraph.Contains(neighbor2)) {
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        for (int const nNeighbor : neighbors[neighbor]) if(vMarkedVertices[nNeighbor]) {
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
        REMOVE_VERTEX(partition, vertex);
        for (const int neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            assert(partitions[neighbor] == partition);
            REMOVE_VERTEX(partition, neighbor);
            remaining.Remove(neighbor);
            independent_set[neighbor] = 1;
            for (const int nNeighbor : neighbors[neighbor]) if(inGraph.Contains(nNeighbor)) {
                // assert(neighbors[nNeighbor].Contains(neighbor));
                REMOVE_NEIGHBOR(partition, nNeighbor, neighbor);
                INSERT_REMAINING(partition, remaining, nNeighbor);
            }
            neighbors[neighbor].Clear();
            assert(!inGraph.Contains(neighbor));
        }
        isolatedCliqueCount += deg + 1;
        neighbors[vertex].Clear();

        // vReductions.emplace_back(std::move(reduction));

        profilingAddTimeSuccessfulIsolatedClique(&profilingHelper, partition);
        return true;
    }
    assert(false);
    return false;
}

bool parallel_reductions::isTwoNeighborhoodInSamePartition(int const vertex, int const partition, ArraySet &remaining) {
    if(isBoundaryVertex(vertex)) {
        return false;
    }
    for(int neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        if(isBoundaryVertex(neighbor)) {
            return false;
        } else if(neighborhoodChanged.Contains(neighbor)) {
            neighborhoodChanged.Remove(neighbor);
            remaining.Insert(neighbor);
        }
    }
    return true;
}

bool parallel_reductions::removeUnconfined(int const partition, int const vertex, ArraySet &remaining, fast_set &closedNeighborhood, vector<int> &neighborhood, vector<int> &numNeighborsInS, vector<int> &neighborsInS, int &removedUnconfinedVerticesCount, int &numDiamondReductions) {
    assert(neighborhood.size() >= neighbors.size());
    assert(numNeighborsInS.size() >= neighbors.size());
    assert(neighborsInS.size() >= 2 * neighbors.size());
    closedNeighborhood.clear();
    closedNeighborhood.add(vertex);
    int sizeS = 1, sizeNeighborhood = 0;
    for (int u : neighbors[vertex]) if(inGraph.Contains(u)) {
        closedNeighborhood.add(u);
        if(partitions[u] == partition) {
            neighborhood[sizeNeighborhood++] = u;
            numNeighborsInS[u] = 1;
        }
    }
    bool vertexAddedToS = true;

    while (vertexAddedToS) {
        vertexAddedToS = false;
        for (int i = 0; i < sizeNeighborhood; i++) {
            int const u = neighborhood[i];
            if (numNeighborsInS[u] != 1)  {
                continue;
            } 
            int neighborToAdd = -1;
            for (int const w : neighbors[u]) if(inGraph.Contains(w) && !closedNeighborhood.get(w)) {
                if (neighborToAdd >= 0) {
                    // There is more than 1 neighbor outside of N[S]
                    neighborToAdd = -2;
                    break;
                }
                neighborToAdd = w;
            }
            if (neighborToAdd == -1) {
                // There is a vertex in N(u) that doesn't have any neighbors outside of N[S]
                // Input vertex is unconfined
                independent_set[vertex] = 1;
                inGraph.Remove(vertex);
                for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                    REMOVE_NEIGHBOR(partition, neighbor, vertex);
                    INSERT_REMAINING(partition, remaining, neighbor);
                }
                neighbors[vertex].Clear();
                remaining.Remove(vertex);
                ++removedUnconfinedVerticesCount;
                return true;
            } else if (neighborToAdd >= 0) {
                // if(sizeS >= MAX_SIZE_UNCONFINED) {
                //     continue;
                // }
                // there is a vertex in N(u) that has exactly one neighbor outside of N[S]
                // that vertex has to be added to S
                if(partitions[neighborToAdd] == partition) {
                    vertexAddedToS = true;
                    closedNeighborhood.add(neighborToAdd);
                    sizeS++;
                    for (int w : neighbors[neighborToAdd]) if(inGraph.Contains(w)) {
                        if (closedNeighborhood.add(w)) {
                            if(partitions[w] == partition) {
                                neighborhood[sizeNeighborhood++] = w;
                                numNeighborsInS[w] = 1;
                            }
                        } else {
                            if(partitions[w] == partition)
                                numNeighborsInS[w]++;
                        }
                    }
                }
            }
        }
    }
    
    if (sizeS >= 2) {
        closedNeighborhood.clear();
        int N = neighbors.size();
        for (int i = 0; i < sizeNeighborhood; i++) closedNeighborhood.add(neighborhood[i]);
        for (int i = 0; i < sizeNeighborhood; i++) {
            neighborsInS[i] = neighborsInS[N + i] = -1;
            int u = neighborhood[i];
            if (numNeighborsInS[u] != 2) continue;
            int v1 = -1, v2 = -1;
            // numNeighborsInS[u] == 2 assures that there are exactly two neighbors in S
            // !closedNeighborhood.get(w) can only cause the loop to find more vertices, not less
            // => only vertices with exactly 2 neighbors outside of N(S) are found
            for (int w : neighbors[u]) if (inGraph.Contains(w) && !closedNeighborhood.get(w)) {
                if (v1 < 0) v1 = w;
                else if (v2 < 0) v2 = w;
                else {
                    v1 = v2 = -1;
                    break;
                }
            }
            if (v1 > v2) {
                int t = v1;
                v1 = v2;
                v2 = t;
            }
            neighborsInS[i] = v1;
            neighborsInS[N + i] = v2;
        }
        for (int i = 0; i < sizeNeighborhood; i++) if (neighborsInS[i] >= 0 && neighborsInS[N + i] >= 0) {
            int u = neighborhood[i];
            closedNeighborhood.clear();
            // TODO
            for (int w : neighbors[u]) if(inGraph.Contains(w)) closedNeighborhood.add(w);
            for (int j = i + 1; j < sizeNeighborhood; j++) if (neighborsInS[i] == neighborsInS[j] && neighborsInS[N + i] == neighborsInS[N + j] && !closedNeighborhood.get(neighborhood[j])) {
                // Vertex is unconfined
                independent_set[vertex] = 1;
                inGraph.Remove(vertex);
                for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                    REMOVE_NEIGHBOR(partition, neighbor, vertex);
                    INSERT_REMAINING(partition, remaining, neighbor);              
                }
                neighbors[vertex].Clear();
                remaining.Remove(vertex);
                ++numDiamondReductions;
                return true;
            }
        }
    }
    return false;
}

bool parallel_reductions::removeTwin(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount)
{
    assert(partitions[vertex] == partition);
    assert(vMarkedVertices.size() == neighbors.size());
    // This takes really long (it's O(n))
    // assert(std::accumulate(vMarkedVertices.begin(), vMarkedVertices.end(), false, std::logical_or<bool>()) == false);
    if(isBoundaryVertex(vertex))
        return false;

    if(degree(vertex) != 3)
        return false;

    int twinNeighbors[3];
    int numNeighbors = 0;


    for(int const neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        assert(numNeighbors < 3);
        twinNeighbors[numNeighbors++] = neighbor;
    }


    int smallestDegreeNeighbor = -1;
    int smallesDegreeNeighborDegree = INT_MAX;
    for(int i = 0; i < 3; ++i) {
        int const neighbor = twinNeighbors[i];
        assert(partitions[neighbor] == partition);
        assert(neighbor != vertex);

        vMarkedVertices[neighbor] = true;
        int const neighborDegree = degree(neighbor);
        if(neighborDegree < smallesDegreeNeighborDegree) {
            smallesDegreeNeighborDegree = neighborDegree;
            smallestDegreeNeighbor = neighbor;
        }
    }
    assert(smallestDegreeNeighbor != -1);

    int twin = -1;
    for(int possibleTwin: neighbors[smallestDegreeNeighbor]) if(inGraph.Contains(possibleTwin)) {
        if(possibleTwin == vertex) continue;
        if(partitions[possibleTwin] != partition) continue;
        if(vMarkedVertices[possibleTwin]) continue;
        if(degree(possibleTwin) != 3) continue;
        assert(partitions[possibleTwin] == partitions[vertex]);
        bool isTwin = true;
        int neighborCount(0);
        for(int twinNeighbor: neighbors[possibleTwin]) if(inGraph.Contains(twinNeighbor)) {
            if(!vMarkedVertices[twinNeighbor]) {
                isTwin = false;
                break;
            }
            ++neighborCount;
        }
        if(isTwin && neighborCount == 3) {
            twin = possibleTwin;
            break;
        }
    }
    if(twin == -1) {
        for(int i = 0; i < 3; ++i) {
            int const markedVertex = twinNeighbors[i];
            vMarkedVertices[markedVertex] = false;
        }
        return false;
    }
    assert(twin >= 0);
    assert(partitions[twin] == partitions[vertex]);
    assert(!isBoundaryVertex(twin));
    assert(!neighbors[vertex].Contains(twin));

    bool isNeighborhoodAdjacent = false;
    for(int i = 0; i < 3; ++i) {
        int const neighbor1 = twinNeighbors[i];
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
        REMOVE_VERTEX(partition, vertex);
        independent_set[vertex] = 0;
        remaining.Remove(vertex);
        REMOVE_VERTEX(partition, twin);
        independent_set[twin] = 0;
        remaining.Remove(twin);
        for(int i = 0; i < 3; ++i) {
            int const neighbor1 = twinNeighbors[i];
            assert(partitions[neighbor1] == partition);
            REMOVE_VERTEX(partition, neighbor1);
            remaining.Remove(neighbor1);
            independent_set[neighbor1] = 1;
            for(int const neighbor2: neighbors[neighbor1]) if(inGraph.Contains(neighbor2)) {
                REMOVE_NEIGHBOR(partition, neighbor2, neighbor1);
                INSERT_REMAINING(partition, remaining, neighbor2);
            }
            neighbors[neighbor1].Clear();
        }
        neighbors[vertex].Clear();
        neighbors[twin].Clear();
        removedTwinCount += 5;
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
            for(int i = 0; i < 3; ++i) {
                int const neighbor = twinNeighbors[i];
                reduction.AddNeighbor(neighbor);
            }

            int neighborHoodSize(0);
            for(int neighbor: reduction.GetNeighbors()) {
                assert(partitions[neighbor] == partitions[twin]);
                // TODO
                neighbors[neighbor].Remove(twin);
                // TODO
                neighbors[neighbor].Remove(vertex);
                neighborHoodSize += degree(neighbor);
            }
            neighbors[twin].Clear();
            neighbors[vertex].Clear();
            neighbors[vertex].Resize(neighborHoodSize);
            for(int neighbor1: reduction.GetNeighbors()) {
                assert(!isBoundaryVertex(neighbor1));
                for(int neighbor2: neighbors[neighbor1]) if(inGraph.Contains(neighbor2)) {
                    assert(partitions[neighbor2] == partitions[neighbor1]);
                    assert(neighbor2 != vertex);
                    assert(neighbor2 != twin);
                    assert(!vMarkedVertices[neighbor2]);
                    // TODO
                    neighbors[neighbor2].Remove(neighbor1);
                    vertexDegree[neighbor2]--;
                    neighbors[vertex].Insert(neighbor2);
                    if(!neighbors[neighbor2].Contains(vertex))
                        vertexDegree[neighbor2]++;
                    neighbors[neighbor2].Insert(vertex);
                    remaining.Insert(neighbor2);
                }
                neighbors[neighbor1].Clear();
                REMOVE_VERTEX(partition, neighbor1);
                remaining.Remove(neighbor1);
            }
            vertexDegree[vertex] = neighbors[vertex].Size();
            REMOVE_VERTEX(partition, twin);
            assert(!isBoundaryVertex(twin));
            remaining.Remove(twin);
            remaining.Insert(vertex);
            vReductions.push_back(reduction);
            assert(inGraph.Contains(vertex));
            assert(!inGraph.Contains(reduction.GetNeighbors()[0]));
            assert(!inGraph.Contains(reduction.GetNeighbors()[1]));
            assert(!inGraph.Contains(reduction.GetNeighbors()[2]));
            foldedTwinCount += 4;
            reduced = true;
        }
    }

    for(int i = 0; i < 3; ++i) {
        int const markedVertex = twinNeighbors[i];
        vMarkedVertices[markedVertex] = false;
    }
    return reduced;
}

bool parallel_reductions::FoldVertex(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, int &foldedVertexCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    if(degree(vertex) != 2) {
        profilingAddTimeUnsuccessfulFoldDegree(&profilingHelper, partition);
        return false;
    }

    if (!isTwoNeighborhoodInSamePartition(vertex, partition, remaining)) { 
        profilingAddTimeUnsuccessfulFoldWrongPartition(&profilingHelper, partition);
        return false;
    }

    int neighbor1 = -1;
    int neighbor2 = -1;
    for(int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        if(neighbor1 == -1)
            neighbor1 = neighbor;
        else if (neighbor2 == -1) 
            neighbor2 = neighbor;
        else {
            assert(false);
        }
    }
    assert(neighbor2 != -1);

    int const vertex1(neighbor1);
    int const vertex2(neighbor2);

    for(int const neighbor2 : neighbors[vertex1]) if(inGraph.Contains(neighbor2)) {
        if(neighbor2 == vertex2) {
            profilingAddTimeUnsuccessfulFoldAdjacent(&profilingHelper, partition);
        return false; // neighbors must not be adjacent.
        }
    }

    foldedVertexCount += 2;

    assert(partitions[vertex1] == partition);
    assert(partitions[vertex2] == partition);

    Reduction reduction(FOLDED_VERTEX);
    reduction.SetVertex(vertex);
    reduction.AddNeighbor(vertex1);
    reduction.AddNeighbor(vertex2);

    neighbors[vertex].Clear();
    int const vertex1degree = degree(vertex1);
    int const vertex2degree = degree(vertex2);
    neighbors[vertex].Resize(vertex1degree + vertex2degree);
    // neighbors[vertex1].Remove(vertex);
    // neighbors[vertex2].Remove(vertex);

    for (int const neighbor1 : neighbors[vertex1]) if(inGraph.Contains(neighbor1)) {
        assert(partitions[neighbor1] == partition);
        if (neighbor1 == vertex) continue;
        // TODO: is it possible to not remove this?
        neighbors[neighbor1].Remove(vertex1);
        vertexDegree[neighbor1]--;
        neighbors[vertex].Insert(neighbor1);
        assert(partitions[vertex] == partitions[neighbor1]);
        INSERT_REMAINING(partition, remaining, neighbor1);
        // remaining.Insert(neighbor1);
    }
    neighbors[vertex1].Clear();
    assert(partitions[vertex] == partitions[vertex1]);

    for (int const neighbor2 : neighbors[vertex2]) if(inGraph.Contains(neighbor2)) {
        assert(partitions[neighbor2] == partition);
        if (neighbor2 == vertex) continue;
        // TODO: is it possible to not remove this?
        neighbors[neighbor2].Remove(vertex2);
        vertexDegree[neighbor2]--;
        neighbors[vertex].Insert(neighbor2);
        assert(partitions[vertex] == partitions[neighbor2]);
        INSERT_REMAINING(partition, remaining, neighbor2);
        // remaining.Insert(neighbor2);
    }

    neighbors[vertex2].Clear();
    assert(partitions[vertex] == partitions[vertex2]);
    
    for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        if(!neighbors[neighbor].Contains(vertex))
            vertexDegree[neighbor]++;
        neighbors[neighbor].Insert(vertex);
    }

    vertexDegree[vertex] = neighbors[vertex].Size();

    INSERT_REMAINING(partition, remaining, vertex);
    // remaining.Insert(vertex);

    vReductions.emplace_back(std::move(reduction));

    remaining.Remove(vertex1);
    remaining.Remove(vertex2);
    REMOVE_VERTEX(partition, vertex2);
    REMOVE_VERTEX(partition, vertex1);

    profilingAddTimeSuccessfulFold(&profilingHelper, partition);
    return true;
}


void parallel_reductions::UpdateRemaining(vector<ArraySet> &remainingPerPartition, vector<vector<int>> &bufferPerPartition) {
    int numPartitions = remainingPerPartition.size();
#pragma omp parallel for schedule(dynamic)
    for(int partition = 0; partition < numPartitions; ++partition) {
      auto tid = omp_get_thread_num();
      vector<int> &buffer = bufferPerPartition[tid];
        ArraySet &remaining = remainingPerPartition[partition];
        int numVerticesRemoved = 0;
        for(int vertex: inGraphPerPartition[partition]) {
            assert(partitions[vertex] == partition);
            if(!inGraph.Contains(vertex)) {
                neighborhoodChanged.Remove(vertex);
                remaining.Remove(vertex);
                buffer[numVerticesRemoved++] = vertex;
            }
            else if(neighborhoodChanged.Contains(vertex)) {
                remaining.Insert(vertex);
            }
        }
        for(int i = 0; i < numVerticesRemoved; ++i) {
            int vertex = buffer[i];
            inGraphPerPartition[partition].Remove(vertex);
        }
    }
}

bool parallel_reductions::LPReduction(vector<ArraySet> &remainingPerPartition, vector<vector<int>> &bufferPerPartition, int &numLPReductions) {
    int sizeBefore = inGraph.Size();
    int N = neighbors.size();
    double startTime = omp_get_wtime();
    UpdateRemaining(remainingPerPartition, bufferPerPartition);
    double updateRemainingBeforeTime = omp_get_wtime();
    maximumMatching.LoadGraph(neighbors, inGraph, vertexDegree);
    double loadGraphTime = omp_get_wtime();
    maximumMatching.KarpSipserInit(inGraph);
    double initTime = omp_get_wtime();
    maximumMatching.MS_BFS_Graft();
    double maximumMatchingTime = omp_get_wtime();
    maximumMatching.MarkReachableVertices();
    double markVerticesTime = omp_get_wtime();
    bool changed = false;
#pragma omp parallel for
    for(int vertex = 0; vertex < N; ++vertex) {
        if(!inGraph.Contains(vertex))
            continue;
        if(maximumMatching.reachableVertices[vertex] == 0 && maximumMatching.reachableVertices[vertex + N] > 0) {
            changed = true;
            // vertex is in the vertex cover
            independent_set[vertex] = 1;
            inGraph.Remove(vertex);
            for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                vertexDegree[neighbor]--;
                if(partitions[vertex] != partitions[neighbor])
                    numCutEdges[neighbor]--;
                neighborhoodChanged.Insert(neighbor);
            }
            neighbors[vertex].Clear();
        } else if (maximumMatching.reachableVertices[vertex] > 0 && maximumMatching.reachableVertices[vertex + N] == 0) {
            changed = true;
            // vertex is in the independent set
            // Nothing to to for the neighbors because they get removed too (two vertices on the same edge can't be 0)
            independent_set[vertex] = 0;
            inGraph.Remove(vertex);
            neighbors[vertex].Clear();
        }
        // else: We don't know it
    }
    double applyReductionTime = omp_get_wtime();

    UpdateRemaining(remainingPerPartition, bufferPerPartition);
    double updateRemainingAfterTime = omp_get_wtime();

    int sizeAfter = inGraph.Size();

/*    std::cout << "Time for UpdateRemaining (before reduction): " << updateRemainingBeforeTime - startTime << std::endl;
    std::cout << "Time for loading the graph: " << loadGraphTime - updateRemainingBeforeTime << std::endl;
    std::cout << "Time for KarpSipserInit: " << initTime - loadGraphTime << std::endl;
    std::cout << "Time for MS_BFS_Graft: " << maximumMatchingTime - initTime << std::endl;
    std::cout << "Time for MarkReachableVertices: " << markVerticesTime - maximumMatchingTime << std::endl;
    std::cout << "Time for applying result: " << applyReductionTime - markVerticesTime << std::endl;
    std::cout << "Time for UpdateRemaining (after reduction): " << updateRemainingAfterTime - applyReductionTime << std::endl;
    std::cout << "Total time: " << updateRemainingAfterTime - startTime << std::endl;
    std::cout << "Vertices removed by LP reduction : " << sizeBefore - sizeAfter << std::endl;
*/    numLPReductions += sizeBefore - sizeAfter;
    return changed;
}

/*void parallel_reductions::reduce_graph_parallel() {

  int numPartitions = partition_nodes.size();
  vector<long> degreeCount(numPartitions, 0);

  long numThreads;
    #pragma omp parallel
  {
    numThreads = omp_get_num_threads();
  }
  std::cout << "num threads: " << numThreads << std::endl;
  std::cout << "num partitions: " << numPartitions << std::endl;


  auto neighbors_per_partition = std::vector<std::vector<SparseArraySet>>(numPartitions);                                                                                                                                                                                   
                                                                                                                                                                                                                                                                              
  for(int partition = 0; partition < numPartitions; ++partition) {                                                                                                                                                                                                            
    neighbors_per_partition[partition] = vector<SparseArraySet>(neighbors.size());                                                                                                                                                                                            
    for(int vertex = 0; vertex < neighbors.size(); ++vertex) {                                                                                                                                                                                                                
      neighbors_per_partition[partition][vertex] = SparseArraySet(neighbors[vertex].Size());                                                                                                                                                                                  
      for(int const neighbor: neighbors[vertex]) {                                                                                                                                                                                                                            
        neighbors_per_partition[partition][vertex].Insert(neighbor);                                                                                                                                                                                                          
      }                                                                                                                                                                                                                                                                       
    }                                                                                                                                                                                                                                                                         
    }

  double startClock = omp_get_wtime();
    #pragma omp parallel for
  for(int partition = 0; partition < numPartitions; ++partition) {
    long degreeCountLocal = 0;
    for(int i = 0; i < 50; ++i) {
      for (int const vertex : partition_nodes[partition]) {
        if(inGraph.Contains(vertex)) {
          for(int const neighbor : neighbors[vertex]) {
            if(inGraph.Contains(neighbor)) {
                degreeCountLocal++;
                }
          }
        }
      }
    }
    std::cout << partition << ": " << degreeCountLocal << std::endl;
  }
  double endClock = omp_get_wtime();

  std::cout << "Time: " << endClock - startClock << std::endl;
  long sum = 0;
  for(int partition = 0; partition < numPartitions; ++partition) {
    std::cout << partition << ": " << degreeCount[partition] << std::endl;
    sum += degreeCount[partition];
  }
  std::cout << "Degree sum: " << sum << std::endl;


}*/

void parallel_reductions::initGlobalBurstEstimator() {
    global_burst_timer = omp_get_wtime();
    global_burst_estimation = -1.0;
    lastFinishedThread = -1;
    terminationFlag = false;;
    firstFinished = -1;
}
double parallel_reductions::finishThreadAndGetEstimatedBurstLength(int tid) {
    global_burst_timer_mutex.lock();
    double current_time = omp_get_wtime();
    bool firstFinisher = lastFinishedThread == -1;
    lastFinishedThread = tid;
    double lastBurstLength = current_time - global_burst_timer;
    if (!firstFinisher) {
        global_burst_estimation = global_burst_estimation <= 0.0 ? lastBurstLength : GLOBAL_BURST_ESTIMATION_ALPHA * lastBurstLength + (1 - GLOBAL_BURST_ESTIMATION_ALPHA) * global_burst_estimation;
    } else {
        firstFinished = tid;
    }
    global_burst_timer = current_time;
    global_burst_timer_mutex.unlock();
    if(firstFinisher) {
        return 0.0;
    } else {
        return global_burst_estimation;
    }
}
bool parallel_reductions::isLastFinishedThread(int tid) {
    return lastFinishedThread == tid && tid != firstFinished;
}
bool parallel_reductions::shouldTerminate() {
    return terminationFlag;
}
void parallel_reductions::terminateOtherThreads() {
    terminationFlag = true;
}

void parallel_reductions::reduce_graph_parallel() {
    long numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    std::cout << "num threads: " << numThreads << std::endl;

    int numPartitions = partition_nodes.size();
    profilingInit(&profilingHelper, &neighbors, numPartitions);

    inGraphPerPartition = vector<ArraySet>(numPartitions);

    vector<vector<bool>> vMarkedVerticesPerTid(numThreads);
    vector<vector<int>> tempInt1PerTid(numThreads);
    vector<vector<int>> tempInt2PerTid(numThreads);
    vector<fast_set> fastSetPerTid(numThreads, fast_set(0));
    vector<vector<int>> tempIntDoubleSizePerTid(numThreads);
    vector<ArraySet> remainingPerPartition(numPartitions);
    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(numPartitions);
    for(int partition = 0; partition < numPartitions; partition++) {
      ArraySet remaining = ArraySet(m_AdjacencyArray.size());
      remainingPerPartition[partition] = remaining;
      inGraphPerPartition[partition] = ArraySet(m_AdjacencyArray.size());
      for (int const vertex : partition_nodes[partition]) {
        if(inGraph.Contains(vertex)) {
          assert(partitions[vertex] == partition);
          inGraphPerPartition[partition].Insert(vertex);
        }
      };
    }
#pragma omp parallel for
    for(int tid = 0; tid < numThreads; tid++) {
        // std::cout << "Start allocating memory for block " << tid << std::endl;
        vMarkedVerticesPerTid[tid] = std::vector<bool>(m_AdjacencyArray.size(), false);
        tempInt1PerTid[tid] = vector<int>(m_AdjacencyArray.size());
        tempInt2PerTid[tid] = vector<int>(m_AdjacencyArray.size());
        fastSetPerTid[tid] = fast_set(m_AdjacencyArray.size());
        tempIntDoubleSizePerTid[tid] = vector<int>(m_AdjacencyArray.size() * 2);
    }
    std::cout << "Finished allocating memory" << std::endl;

    vector<double> partitionTimes(numPartitions);
    vector<double> partitionFinishTimes(numPartitions);
    vector<int> partitionFinishSizes(numPartitions);
    vector<int> numIsolatedCliqueReductions(numPartitions, 0);
    vector<int> numVertexFoldReductions(numPartitions, 0);
    vector<int> numTwinReductionsRemoved(numPartitions, 0);
    vector<int> numTwinReductionsFolded(numPartitions, 0);
    vector<int> removedUnconfinedVerticesCount(numPartitions, 0);
    vector<int> numDiamondReductions(numPartitions, 0);
    dependecy_checking_burst_estimation = std::vector<double>(numPartitions, -1.0);
    dependency_checking_times = std::vector<double>(numPartitions, 0.0);

    global_burst_timer = 0.0;
    global_burst_estimation = -1.0;
    lastFinishedThread = -1;
    terminationFlag = false;

    int numLPReductions = 0;
    
    assert(checkDegrees());

    double startClock = omp_get_wtime();
    double tmpClock;
    double LPTime = 0;
    double restTime = 0;
    std::cout << "Filling remaining vertices" << std::endl;
#pragma omp parallel for schedule(dynamic)
    for(int partition = 0; partition < numPartitions; ++partition) {
        remainingPerPartition[partition].Clear();
        for (int const vertex : partition_nodes[partition]) {
            if(inGraph.Contains(vertex)) {
                assert(partitions[vertex] == partition);
                remainingPerPartition[partition].Insert(vertex);
            }
        }
    }


    // std::cout << "Start LP reduction" << std::endl;
    // tmpClock = omp_get_wtime();

    // LPReduction(remainingPerPartition, tempInt1PerTid, numLPReductions);
    // LPTime += omp_get_wtime() - tmpClock;

    // std::cout << "done with LP reduction" << std::endl;

    bool changed = true;
    int numIterations = 0;
    double terminationTime = 0.0;
    while(changed) {
        terminationTime = -1.0;
        initGlobalBurstEstimator();
        std::cout << "Iteration " << numIterations << " starts at " << omp_get_wtime() - startClock << " current size: " << inGraph.Size() << std::endl;
        // int sizeBefore = inGraph.Size();
      //   std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
      // std::cout << "starting new iteration: " << numIterations << std::endl;
      // std::cout << "Current rest time: " << restTime << std::endl;
      // auto partitionTimesCopy(partitionTimes);
      // auto old_graphsize = inGraph.Size();
        int graphSizeLastSample = 0;
        for(int partition = 0; partition < numPartitions; ++partition) {
            graphSizeLastSample += inGraphPerPartition[partition].Size();
        }
        double timeLastSample = omp_get_wtime();
        double timeBefore = timeLastSample;
        int graphSizeBefore = graphSizeLastSample;

        tmpClock = omp_get_wtime();
#pragma omp parallel for schedule(dynamic,1)
        for(int partition = 0; partition < numPartitions; partition++) {
          auto tid = omp_get_thread_num();
          // std::cout << "partition " << partition << " on tid " << tid << std::endl;
          //std::cout << partition << ": starting new iteration" << std::endl;
            ApplyReductions(partition, ReductionsPerPartition[partition], vMarkedVerticesPerTid[tid], remainingPerPartition[partition], tempInt1PerTid[tid], tempInt2PerTid[tid], fastSetPerTid[tid], tempIntDoubleSizePerTid[tid], partitionTimes[partition], numIsolatedCliqueReductions[partition], numVertexFoldReductions[partition], numTwinReductionsRemoved[partition], numTwinReductionsFolded[partition], removedUnconfinedVerticesCount[partition], numDiamondReductions[partition]);
            partitionFinishTimes[partition] = omp_get_wtime() - startClock;
            partitionFinishSizes[partition] = inGraph.Size();
            // if(!shouldTerminate()) {
            //     double sleeptime = finishThreadAndGetEstimatedBurstLength(tid) * GLOBAL_BURST_THRESHOLD_MULTIPLIER;
            //     double startTime = omp_get_wtime();
            //     std::cout << tid << " sleeping for " << sleeptime << std::endl;
            //     while( (omp_get_wtime() - startTime) < sleeptime);
            //     if(isLastFinishedThread(tid)) {
            //         std::cout << tid << " Terminating others" << std::endl;
            //         terminationTime = omp_get_wtime() - startClock;
            //         terminateOtherThreads();
            //     }
            // }

            global_burst_timer_mutex.lock();
            bool isFirstFinisher = firstFinished == -1;
            firstFinished = 1;
            global_burst_timer_mutex.unlock();
            if(isFirstFinisher) {
                int graphSizeCurrentSample = 0;
                for(int partition = 0; partition < numPartitions; ++partition) {
                    graphSizeCurrentSample += inGraphPerPartition[partition].Size();
                }
                double current_time = omp_get_wtime();
                // double last_delta = (graphSizeLastSample - graphSizeCurrentSample) / (current_time - timeLastSample);
                double last_delta = 0;
                graphSizeLastSample = graphSizeCurrentSample;
                timeLastSample = current_time;
                while(true) {
                    double startTime = omp_get_wtime();
                    while( (omp_get_wtime() - startTime) < 0.2);
                    graphSizeCurrentSample = 0;
                    for(int partition = 0; partition < numPartitions; ++partition) {
                        graphSizeCurrentSample += inGraphPerPartition[partition].Size();
                    }
                    current_time = omp_get_wtime();
                    double current_delta = (graphSizeLastSample - graphSizeCurrentSample) / (current_time - timeLastSample);
                    double global_delta = (graphSizeBefore - graphSizeCurrentSample) / (current_time - timeBefore);
                    graphSizeLastSample = graphSizeCurrentSample;
                    timeLastSample = current_time;
                    if(current_delta <= 0.001 * global_delta) {
                        terminationTime = omp_get_wtime() - startClock;
                        terminateOtherThreads();
                        break;
                    }
                    last_delta = current_delta;
                }
            }
        }
        restTime += omp_get_wtime() - tmpClock;
        if(terminationTime > 0.0) {
            std::cout << "Termination time: " << terminationTime << std::endl;
        }
        for(int partition = 0; partition < numPartitions; partition++) {
            std::cout << "Partition " << partition << " finished iteration " << numIterations << " at " << partitionFinishTimes[partition] << " with (approximate) size " << partitionFinishSizes[partition] << std::endl;
        }
        // int sizeAfter = inGraph.Size();
        // std::cout << "Vertices removed by other reductions: " << sizeBefore - sizeAfter << std::endl;
        // changed = false;
        tmpClock = omp_get_wtime();

        changed = LPReduction(remainingPerPartition, tempInt1PerTid, numLPReductions);
        // changed = false;
        LPTime += omp_get_wtime() - tmpClock;
        // std::cout << "Size after iteration: " << inGraph.Size() << std::endl;
        std::cout << "Graph size after iteration " << numIterations << " at time " << omp_get_wtime() - startClock << ": " << inGraph.Size();
        std::cout << " -- Current rest time: " << restTime << " -- Current LP time: " << LPTime << std::endl;
        numIterations++;
    }
    std::cout << "Num iterations: " << numIterations << std::endl;


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

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Number of unconfined vertices removed: " << removedUnconfinedVerticesCount[partition] << endl;
    }

    for(int partition = 0; partition< numPartitions; partition++) {
        cout << partition << ": Number of diamond reductions: " << numDiamondReductions[partition] << endl;
    }

    cout << "Number of vertices removed by LP reduction: " << numLPReductions << endl;
    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;

    std::cout << "LP time: " << LPTime << std::endl;
    std::cout << "Rest time: " << restTime << std::endl;

    int sum_isolated_clique = std::accumulate(numIsolatedCliqueReductions.begin(), numIsolatedCliqueReductions.end(), 0);
    int sum_vertex_fold = std::accumulate(numVertexFoldReductions.begin(), numVertexFoldReductions.end(), 0);
    int sum_twin_removed = std::accumulate(numTwinReductionsRemoved.begin(), numTwinReductionsRemoved.end(), 0);
    int sum_twin_folded = std::accumulate(numTwinReductionsFolded.begin(), numTwinReductionsFolded.end(), 0);
    int sum_unconfined = std::accumulate(removedUnconfinedVerticesCount.begin(), removedUnconfinedVerticesCount.end(), 0);
    int sum_diamond = std::accumulate(numDiamondReductions.begin(), numDiamondReductions.end(), 0);
    int sum_reductions = sum_isolated_clique + sum_vertex_fold + sum_twin_removed + sum_twin_folded + sum_unconfined + sum_diamond + numLPReductions;
    assert(sum_reductions == neighbors.size() - inGraph.Size());
    assert(checkDegrees());
}

bool parallel_reductions::checkDegrees() {
    for(int vertex = 0; vertex < neighbors.size(); ++vertex) if(inGraph.Contains(vertex)) {
        int deg = 0;
        int cutEdges = 0;
        for(int const neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            deg++;
            if(partitions[vertex] != partitions[neighbor])
                cutEdges++;
        }
        if(vertexDegree[vertex] != deg) {
            std::cout << "Vertex " << vertex << " has degree " << deg << " but vertexDegree[vertex] has value " << vertexDegree[vertex] << " number of elements in neighbors[vertex]: " << neighbors[vertex].Size() << std::endl;
            return false;
        }
        if(numCutEdges[vertex] != cutEdges){
            std::cout << "Vertex " << vertex << " has " << cutEdges << " cut edges but numCutEdges[vertex] has value " << numCutEdges[vertex] << std::endl;
            return false;
        }
    }
    return true;
}

void parallel_reductions::reduce_graph_sequential() {
    return;
    long numThreads;

    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    //omp_set_num_threads(1);
    std::cout << "numThreads: " << numThreads << std::endl;
    profilingInit(&profilingHelper, &neighbors, 1);

    int N = m_AdjacencyArray.size();

    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(1);
    std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
    vector<int> tempInt1 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempInt2 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempIntDoubleSize = vector<int>(m_AdjacencyArray.size() * 2);
    fast_set fastSet(m_AdjacencyArray.size());
    ArraySet remaining = ArraySet(N);

    partition_nodes.resize(1);
    partition_nodes[0].clear();
    inGraphPerPartition = vector<ArraySet>(1);
    inGraphPerPartition[0] = ArraySet(N);
    for(int vertex = 0; vertex < N; ++vertex) {
        numCutEdges[vertex] = 0;
        partitions[vertex] = 0;
        partition_nodes[0].push_back(vertex);
        if(inGraph.Contains(vertex)) {
            inGraphPerPartition[0].Insert(vertex);
        }
    }

    double time(0);
    int numIsolatedCliqueReductions(0);
    int numVertexFoldReductions(0);
    int numTwinReductionsRemoved(0);
    int numTwinReductionsFolded(0);
    int removedUnconfinedVerticesCount(0);
    int numLPReductions(0);
    int numDiamondReductions(0);
    
    double startClock = omp_get_wtime();

    vector<ArraySet> remainingPerPartition = {remaining};
    vector<vector<int>> tempInt1PerPartition = {tempInt1};

    std::cout << "Filling remaining vertices set..." << std::endl;
    remaining.Clear();
    for(int vertex = 0; vertex < N; ++vertex)  {
        if(inGraph.Contains(vertex)) {
            assert(partitions[vertex] == 0);
            remaining.Insert(vertex);
        }
    }

    bool changed = true;
    int numIterations = 0;
    while(changed) {
        numIterations++;
        ApplyReductions(0, ReductionsPerPartition[0], vMarkedVertices, remaining, tempInt1, tempInt2, fastSet, tempIntDoubleSize, time, numIsolatedCliqueReductions, numVertexFoldReductions, numTwinReductionsRemoved, numTwinReductionsFolded, removedUnconfinedVerticesCount, numDiamondReductions);
        changed = LPReduction(remainingPerPartition, tempInt1PerPartition, numLPReductions);
    }

    double endClock = omp_get_wtime();
    std::cout << "Num iterations: " << numIterations << std::endl;
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
    cout << "Number of isolated clique reductions: " << numIsolatedCliqueReductions << endl;
    cout << "Number of vertex fold reductions: " << numVertexFoldReductions << endl;
    cout << "Number of twin reductions (removed): " << numTwinReductionsRemoved << endl;
    cout << "Number of twin reductions (folded): " << numTwinReductionsFolded << endl;
    cout << "Number of unconfined vertices removed: " << removedUnconfinedVerticesCount << endl;
    cout << "Number of diamond reductions: " << numDiamondReductions << endl;
    cout << "Number of vertices removed by LP reduction: " << numLPReductions << endl;
    omp_set_num_threads(numThreads);
}

void parallel_reductions::initDependencyCheckingEstimation(int partition) {
    dependency_checking_times[partition] = omp_get_wtime();
}

void parallel_reductions::updateDependencyCheckingEstimation(int partition) {
    double current_time = omp_get_wtime();
    double timeSinceLastReduction = current_time - dependency_checking_times[partition];
    dependecy_checking_burst_estimation[partition] = dependecy_checking_burst_estimation[partition] <= 0.0 ? timeSinceLastReduction : DEPENDENCYCHECKING_BURST_ESTIMATION_ALPHA * timeSinceLastReduction + (1 - DEPENDENCYCHECKING_BURST_ESTIMATION_ALPHA) * dependecy_checking_burst_estimation[partition];
    dependency_checking_times[partition] = current_time;
}

bool parallel_reductions::shouldStopDependencyCheckingReductions(int partition) {
    return dependecy_checking_burst_estimation[partition] <= 0.0 ? false : ((omp_get_wtime() - dependency_checking_times[partition]) > dependecy_checking_burst_estimation[partition] * DEPENDENCY_CHECKING_THRESHOLD_MULTIPLIER);
}

void parallel_reductions::ApplyReductions(int const partition, vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining, vector<int> &tempInt1, vector<int> &tempInt2, fast_set &fastSet, vector<int> &tempIntDoubleSize, double &time, int &isolatedCliqueCount, int &foldedVertexCount, int &removedTwinCount, int &foldedTwinCount, int &removedUnconfinedVerticesCount, int &numDiamondReductions)
{
    double startClock = omp_get_wtime();
    int dependencyCheckingIterations(0);
    int nonDependencyCheckingIterations(0);
    // std::cout << partition << ": Starting reductions..." << std::endl;
    bool changed = true;
    while (changed) {
        changed = false;
      //std::cout << partition << ": Starting reductions with dependency checking..." << std::endl;
        initDependencyCheckingEstimation(partition);
        while (!remaining.Empty()) {
            if(shouldTerminate()) {
                break;
            }
            int const vertex = *(remaining.begin());
            remaining.Remove(vertex);
            assert(inGraph.Contains(vertex));
            assert(partitions[vertex] == partition);
            assert(independent_set[vertex] == -1);
            bool reduction = RemoveIsolatedClique(partition, vertex, vReductions, remaining, vMarkedVertices, isolatedCliqueCount);
            if (!reduction && m_bAllowVertexFolds) {
                reduction = FoldVertex(partition, vertex, vReductions, remaining, foldedVertexCount);
		}
            if (!reduction) {
                reduction = removeTwin(partition, vertex, vReductions, remaining, vMarkedVertices, removedTwinCount, foldedTwinCount);
		}
            // if (!reduction) {
            //     reduction = removeUnconfined(partition, vertex, remaining, fastSet, tempInt1, tempInt2, tempIntDoubleSize, removedUnconfinedVerticesCount, numDiamondReductions);
            //     if(reduction) {
            //         inGraphPerPartition[partition].Remove(vertex);
            //     }
            // }
            /*dependencyCheckingIterations++;
            if(dependencyCheckingIterations % 1000000 == 0) {
                std::cout << partition << ": " << dependencyCheckingIterations << " iterations. Currently queued vertices: " << remaining.Size() << ". Isolated clique reductions: " << isolatedCliqueCount << ", vertex fold count: " << foldedVertexCount << ", twin reduction count (removed): " << removedTwinCount  << ", twin reduction count (folded): " << foldedTwinCount << std::endl;
		}*/
            // if(reduction) {
            //     // Update burst estimation
            //     updateDependencyCheckingEstimation(partition);
            // } else {
            //     // Check if last reduction was long ago and stop dependency checking reductions if threshold is exceeded
            //     if(shouldStopDependencyCheckingReductions(partition)) {
            //         break;
            //     }
            // }
        }
        // std::cout << partition << ": " << dependencyCheckingIterations << " iterations. Isolated clique reductions: " << isolatedCliqueCount << ", vertex fold count: " << foldedVertexCount << ", twin reduction count (removed): " << removedTwinCount  << ", twin reduction count (folded): " << foldedTwinCount << std::endl;
        // std::cout << partition << ": Starting reductions without dependency checking..." << std::endl;
        std::vector<int> verticesToRemove;
        for (int const vertex : inGraphPerPartition[partition]) {
            if(shouldTerminate()) {
                break;
            }
            assert(partitions[vertex] == partition);
            if(inGraph.Contains(vertex)) {
                bool reduction = removeUnconfined(partition, vertex, remaining, fastSet, tempInt1, tempInt2, tempIntDoubleSize, removedUnconfinedVerticesCount, numDiamondReductions);
                if(reduction) {
                    changed = true;
                    verticesToRemove.push_back(vertex);
                }
            }
            nonDependencyCheckingIterations++;
        }
        for(int vertex : verticesToRemove) {
            inGraphPerPartition[partition].Remove(vertex);
        }
        // remaining.Clear();
        // std::cout << partition << ": " << nonDependencyCheckingIterations << " iterations. Unconfined reductions: " << removedUnconfinedVerticesCount << std::endl;

        if(shouldTerminate()) {
            changed = false;
        }
	}
    // std::cout << partition << ": Finished reductions!" << std::endl;
    double endClock = omp_get_wtime();
    time += (endClock - startClock);
}

void parallel_reductions::reduce_graph_sequential_reduction_wise() {
    long numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    omp_set_num_threads(1);
    profilingInit(&profilingHelper, &neighbors, 1);

    int N = m_AdjacencyArray.size();

    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(1);
    std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
    vector<int> tempInt1 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempInt2 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempIntDoubleSize = vector<int>(m_AdjacencyArray.size() * 2);
    fast_set fastSet(m_AdjacencyArray.size());
    ArraySet remainingUse = ArraySet(N);
    ArraySet remaining2 = ArraySet(N);

    partition_nodes.resize(1);
    partition_nodes[0].clear();
    inGraphPerPartition = vector<ArraySet>(1);
    inGraphPerPartition[0] = ArraySet(N);
    for(int vertex = 0; vertex < N; ++vertex) {
        numCutEdges[vertex] = 0;
        partitions[vertex] = 0;
        partition_nodes[0].push_back(vertex);
        if(inGraph.Contains(vertex)) {
            inGraphPerPartition[0].Insert(vertex);
        }
    }

    double time(0);
    int numIsolatedCliqueReductions(0);
    int numVertexFoldReductions(0);
    int numTwinReductionsRemoved(0);
    int numTwinReductionsFolded(0);
    int removedUnconfinedVerticesCount(0);
    int numLPReductions(0);
    int numDiamondReductions(0);
    
    double startClock = omp_get_wtime();

    vector<vector<int>> tempInt1PerPartition = {tempInt1};

    remainingUse.Clear();
    for(int vertex = 0; vertex < N; ++vertex)  {
        if(inGraph.Contains(vertex)) {
            assert(partitions[vertex] == 0);
            remainingUse.Insert(vertex);
        }
    }

    ArraySet *remainingUseptr = &remainingUse;
    ArraySet *remainingInsertptr = &remaining2;
    double isolated_clique_time = 0.0;
    double unconfined_time = 0.0;
    double lp_time = 0.0;
    double vertex_fold_time = 0.0;
    double twin_time = 0.0;
    while(true) {
        bool changed = false;
        ArraySet *temp;
        double start_time;

        start_time = omp_get_wtime();
        changed = RemoveAllIsolatedClique(0, ReductionsPerPartition[0], remainingUseptr, remainingInsertptr, vMarkedVertices, numIsolatedCliqueReductions);
        isolated_clique_time += omp_get_wtime() - start_time;
        temp = remainingInsertptr;
        remainingInsertptr = remainingUseptr;
        remainingUseptr = temp;
        // if(changed) continue;

        start_time = omp_get_wtime();
        changed = removeAllUnconfined(0, remainingUseptr, fastSet, tempInt1, tempInt2, tempIntDoubleSize, removedUnconfinedVerticesCount, numDiamondReductions);
        unconfined_time += omp_get_wtime() - start_time;
        if(changed) continue;

        vector<ArraySet> remainingPerPartition = {*remainingUseptr};
        start_time = omp_get_wtime();
        changed = LPReduction(remainingPerPartition, tempInt1PerPartition, numLPReductions);
        lp_time += omp_get_wtime() - start_time;
        if(changed) continue;

        start_time = omp_get_wtime();
        changed = FoldAllVertices(0, ReductionsPerPartition[0], remainingUseptr, remainingInsertptr, numVertexFoldReductions);
        vertex_fold_time += omp_get_wtime() - start_time;
        temp = remainingInsertptr;
        remainingInsertptr = remainingUseptr;
        remainingUseptr = temp;
        if(changed) continue;

        start_time = omp_get_wtime();
        changed = removeAllTwin(0, ReductionsPerPartition[0], remainingUseptr, remainingInsertptr, vMarkedVertices, numTwinReductionsRemoved, numTwinReductionsFolded);
        twin_time += omp_get_wtime() - start_time;
        temp = remainingInsertptr;
        remainingInsertptr = remainingUseptr;
        remainingUseptr = temp;
        if(changed) continue;

        break;
    }

    double endClock = omp_get_wtime();
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
    cout << "Number of isolated clique reductions: " << numIsolatedCliqueReductions << endl;
    cout << "Number of vertex fold reductions: " << numVertexFoldReductions << endl;
    cout << "Number of twin reductions (removed): " << numTwinReductionsRemoved << endl;
    cout << "Number of twin reductions (folded): " << numTwinReductionsFolded << endl;
    cout << "Number of unconfined vertices removed: " << removedUnconfinedVerticesCount << endl;
    cout << "Number of diamond reductions: " << numDiamondReductions << endl;
    cout << "Number of vertices removed by LP reduction: " << numLPReductions << endl;

    cout << "Isolated clique time: " << isolated_clique_time << std::endl;
    cout << "Unconfined time: " << unconfined_time << std::endl;
    cout << "LP time: " << lp_time << std::endl;
    cout << "Vertex fold time: " << vertex_fold_time << std::endl;
    cout << "Twin time: " << twin_time << std::endl;
    omp_set_num_threads(numThreads);
}

bool parallel_reductions::removeAllUnconfined(int const partition, ArraySet *remainingInsert, fast_set &closedNeighborhood, std::vector<int> &neighborhood, std::vector<int> &numNeighborsInS, std::vector<int> &neighborsInS, int &removedUnconfinedVerticesCount, int &numDiamondReductions) {
    std::vector<int> verticesToRemove;
    bool reduced = false;
    for (int const vertex : inGraphPerPartition[partition]) {
        assert(partitions[vertex] == partition);
        if(inGraph.Contains(vertex)) {
            bool reduction = removeUnconfined(partition, vertex, *remainingInsert, closedNeighborhood, neighborhood, numNeighborsInS, neighborsInS, removedUnconfinedVerticesCount, numDiamondReductions);
            if(reduction) {
                reduced = true;
                verticesToRemove.push_back(vertex);
            }
        }
    }
    for(int vertex : verticesToRemove) {
        inGraphPerPartition[partition].Remove(vertex);
    }
    return reduced;
}
bool parallel_reductions::removeAllTwin(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, std::vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount) {
    bool reduced = false;
    while (!(*remainingUse).Empty()) {
        int const vertex = *((*remainingUse).begin());
        (*remainingUse).Remove(vertex);
        if(!inGraph.Contains(vertex))
            continue;
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        bool reduction = removeTwin(partition, vertex, vReductions, *remainingInsert, vMarkedVertices, removedTwinCount, foldedTwinCount);
        if(!reduction) {
            remainingInsert->Insert(vertex);
        }
        reduced |=reduction;
    }
    return reduced;
}

bool parallel_reductions::RemoveAllIsolatedClique(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, std::vector<bool> &vMarkedVertices, int &isolatedCliqueCount){
    bool reduced = false;
    while (!(*remainingUse).Empty()) {
        int const vertex = *((*remainingUse).begin());
        (*remainingUse).Remove(vertex);
        if(!inGraph.Contains(vertex))
            continue;
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        bool reduction = RemoveIsolatedClique(partition, vertex, vReductions, *remainingInsert, vMarkedVertices, isolatedCliqueCount);
        if(!reduction) {
            remainingInsert->Insert(vertex);
        }
        reduced |=reduction;
    }
    return reduced;
}

bool parallel_reductions::FoldAllVertices(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, int &foldedVertexCount) {
    bool reduced = false;
    while (!(*remainingUse).Empty()) {
        int const vertex = *((*remainingUse).begin());
        (*remainingUse).Remove(vertex);
        if(!inGraph.Contains(vertex))
            continue;
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        bool reduction = FoldVertex(partition, vertex, vReductions, *remainingInsert, foldedVertexCount);
        if(!reduction) {
            remainingInsert->Insert(vertex);
        }
        reduced |=reduction;
    }
    return reduced;
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
