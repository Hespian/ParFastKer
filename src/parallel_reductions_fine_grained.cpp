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

#include "parallel_reductions_fine_grained.h"
#include "ArraySet.h"
#include "SparseArraySet.h"

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
#include <string.h>
#include <memory>

#define ISOLATED_CLIQUE_MAX_NEIGHBORS 2
#define maxTempSize 8192

using namespace std;


parallel_reductions_fine_grained::parallel_reductions_fine_grained(vector<vector<int>> const &adjacencyArray)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size(), true)
 , independent_set(adjacencyArray.size(), -1)
 , maximumMatching(adjacencyArray)
 , vertexDegree(adjacencyArray.size())
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
    for(int node = 0; node < N; ++node) {
        vertexDegree[node] = neighbors[node].Size();
    }
    std::cout << "Finished constructor" << std::endl;
}

parallel_reductions_fine_grained::~parallel_reductions_fine_grained()
{

#ifdef TIMERS
    cout << "Total time spent undoing  reductions  : " << (replaceTimer/(double)CLOCKS_PER_SEC) << endl;
#endif // TIMERS
}

std::vector<std::vector<int>> parallel_reductions_fine_grained::getKernel() {
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

void parallel_reductions_fine_grained::applyKernelSolution(std::vector<int> kernel_solution){
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

int parallel_reductions_fine_grained::degree(int const vertex) {
    return vertexDegree[vertex];
}


bool parallel_reductions_fine_grained::RemoveIsolatedClique(vector<vector<bool>> &vMarkedVerticesPerThread, vector<int> &toRemove, vector<unique_ptr<vector<int>>> &temp)
{
    int const N = neighbors.size();
    int toRemoveCount = 0;
    int numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    vector<int> tempSize = vector<int>(numThreads, 0);
    #pragma omp parallel for
    for(int vertex = 0; vertex < N; vertex++) if(inGraph.Contains(vertex)) {
        {
            int tid = omp_get_thread_num();
            vector<bool> &vMarkedVertices = vMarkedVerticesPerThread[tid];
            int const deg = degree(vertex);

            if(deg > ISOLATED_CLIQUE_MAX_NEIGHBORS)
                goto nextVertex;

            for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                if (degree(neighbor) < deg) {
                    goto nextVertex;
                } else if (degree(neighbor) == deg && neighbor < vertex) {
                    // If it is an isolated clique, the other vertex should be reduced
                    goto nextVertex;
                }
            }

            bool superSet(true);

            for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
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
                    goto nextVertex;
                }
            }
            if (superSet) {
                (*(temp[tid])).push_back(vertex);
                if((*(temp[tid])).size() == maxTempSize) {
                    int startindex = __sync_fetch_and_add(&toRemoveCount, maxTempSize);
                    memcpy(&(toRemove[startindex]), &((*(temp[tid]))[0]), maxTempSize * sizeof(int));
                    (*(temp[tid])).clear();
                }
            }
        }
        nextVertex: ;
    }

    #pragma omp parallel for
    for(int tid = 0; tid < numThreads; ++tid) {
        int startindex = __sync_fetch_and_add(&toRemoveCount, (*(temp[tid])).size());
        memcpy(&(toRemove[startindex]), &((*(temp[tid]))[0]), (*(temp[tid])).size() * sizeof(int));
        (*(temp[tid])).clear();
    }
    // std::cout << "Done with finding isolated cliques! ToRemove: " << toRemoveCount << std::endl;


    #pragma omp parallel for
    for(int i = 0; i < toRemoveCount; ++i) {
        int vertex = toRemove[i];
        assert(inGraph.Contains(vertex));
        independent_set[vertex] = 0;
        inGraph.Remove(vertex);
        for (const int neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            inGraph.Remove(neighbor);
            independent_set[neighbor] = 1;
            for (const int nNeighbor : neighbors[neighbor]) if(inGraph.Contains(nNeighbor)) {
                vertexDegree[nNeighbor]--;
            }
        }
    }

    return toRemoveCount > 0;
}



bool parallel_reductions_fine_grained::LPReduction() {
    int N = neighbors.size();
    double startTime = omp_get_wtime();
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
*/  return changed;
}

bool parallel_reductions_fine_grained::RemoveUnconfined(vector<fast_set> &closedNeighborhoodPerThread, vector<vector<int>> &neighborhoodPerThread, vector<vector<int>> &numNeighborsInSPerThread, vector<vector<int>> &neighborsInSPerThread, vector<char> &isCandidate, vector<int> &candidates, vector<int> &toRemove, vector<unique_ptr<vector<int>>> &temp) {
    int candidateCount = 0;
    int toRemoveCount = 0;

    int numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }

    int const N = neighbors.size();
    #pragma omp parallel for
    for(int vertex = 0; vertex < N; ++vertex) if(inGraph.Contains(vertex)) {
        {
            int tid = omp_get_thread_num();
            fast_set &closedNeighborhood = closedNeighborhoodPerThread[tid];
            vector<int> &neighborhood = neighborhoodPerThread[tid];
            vector<int> &numNeighborsInS = numNeighborsInSPerThread[tid];
            vector<int> &neighborsInS = neighborsInSPerThread[tid];

            closedNeighborhood.clear();
            closedNeighborhood.add(vertex);
            int sizeS = 1, sizeNeighborhood = 0;
            for (int u : neighbors[vertex]) if(inGraph.Contains(u)) {
                closedNeighborhood.add(u);
                neighborhood[sizeNeighborhood++] = u;
                numNeighborsInS[u] = 1;
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
                        (*(temp[tid])).push_back(vertex);
                        if((*(temp[tid])).size() == maxTempSize) {
                            int startindex = __sync_fetch_and_add(&candidateCount, maxTempSize);
                            memcpy(&(candidates[startindex]), &((*(temp[tid]))[0]), maxTempSize * sizeof(int));
                            (*(temp[tid])).clear();
                        }
                        isCandidate[vertex] = true;
                        goto nextVertexFirstRun;
                    } else if (neighborToAdd >= 0) {
                        // there is a vertex in N(u) that has exactly one neighbor outside of N[S]
                        // that vertex has to be added to S
                        vertexAddedToS = true;
                        closedNeighborhood.add(neighborToAdd);
                        sizeS++;
                        for (int w : neighbors[neighborToAdd]) if(inGraph.Contains(w)) {
                            if (closedNeighborhood.add(w)) {
                                neighborhood[sizeNeighborhood++] = w;
                                numNeighborsInS[w] = 1;
                            } else {
                                numNeighborsInS[w]++;
                            }
                        }
                        
                    }
                }
            }
        }
    nextVertexFirstRun:     ;
    }
    #pragma omp parallel for
    for(int tid = 0; tid < numThreads; ++tid) {
        int startindex = __sync_fetch_and_add(&candidateCount, (*(temp[tid])).size());
        memcpy(&(candidates[startindex]), &((*(temp[tid]))[0]), (*(temp[tid])).size() * sizeof(int));
        (*(temp[tid])).clear();
    }
    // std::cout << "Done with run 1! Candidates: " << candidateCount << std::endl;


    #pragma omp parallel for
    for(int i = 0; i < candidateCount; ++i) if(inGraph.Contains(candidates[i])) {
        {
            int vertex = candidates[i];
            int tid = omp_get_thread_num();
            fast_set &closedNeighborhood = closedNeighborhoodPerThread[tid];
            vector<int> &neighborhood = neighborhoodPerThread[tid];
            vector<int> &numNeighborsInS = numNeighborsInSPerThread[tid];
            vector<int> &neighborsInS = neighborsInSPerThread[tid];

            closedNeighborhood.clear();
            closedNeighborhood.add(vertex);
            int sizeS = 1, sizeNeighborhood = 0;
            for (int u : neighbors[vertex]) if(inGraph.Contains(u)) {
                closedNeighborhood.add(u);
                neighborhood[sizeNeighborhood++] = u;
                numNeighborsInS[u] = 1;
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
                        if(isCandidate[u] && u < vertex) {
                            continue;
                        }
                        (*(temp[tid])).push_back(vertex);
                        if((*(temp[tid])).size() == maxTempSize) {
                            int startindex = __sync_fetch_and_add(&toRemoveCount, maxTempSize);
                            memcpy(&(toRemove[startindex]), &((*(temp[tid]))[0]), maxTempSize * sizeof(int));
                            (*(temp[tid])).clear();
                        }
                        goto nextVertexSecond;
                    } else if (neighborToAdd >= 0) {
                        // there is a vertex in N(u) that has exactly one neighbor outside of N[S]
                        // that vertex has to be added to S
                        if(isCandidate[neighborToAdd] && neighborToAdd < vertex) {
                            continue;
                        }
                        if(isCandidate[u] && u < vertex) {
                            continue;
                        }
                        vertexAddedToS = true;
                        closedNeighborhood.add(neighborToAdd);
                        sizeS++;
                        for (int w : neighbors[neighborToAdd]) if(inGraph.Contains(w)) {
                            if (closedNeighborhood.add(w)) {
                                neighborhood[sizeNeighborhood++] = w;
                                numNeighborsInS[w] = 1;
                            } else {
                                numNeighborsInS[w]++;
                            }
                        }
                        
                    }
                }
            }
        }
    nextVertexSecond:     ;
    }
    #pragma omp parallel for
    for(int tid = 0; tid < numThreads; ++tid) {
        int startindex = __sync_fetch_and_add(&toRemoveCount, (*(temp[tid])).size());
        memcpy(&(toRemove[startindex]), &((*(temp[tid]))[0]), (*(temp[tid])).size() * sizeof(int));
        (*(temp[tid])).clear();
    }
    // std::cout << "Done with run 2! ToRemove: " << toRemoveCount << std::endl;

    #pragma omp parallel for
    for(int i = 0; i < candidateCount; ++i) {
        int vertex = candidates[i];
        isCandidate[vertex] = false;
    }  

    #pragma omp parallel for
    for(int i = 0; i < toRemoveCount; ++i) {
        int vertex = toRemove[i];
        independent_set[vertex] = 1;
        inGraph.Remove(vertex);
        for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            vertexDegree[neighbor]--;
        }
        neighbors[vertex].Clear();
    }   
    return toRemoveCount > 0;     
}

void parallel_reductions_fine_grained::reduce_graph_parallel() {
    long numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    std::cout << "num threads: " << numThreads << std::endl;


    vector<vector<int>> tempInt1PerThread(numThreads);
    vector<vector<int>> tempInt2PerThread(numThreads);
    vector<fast_set> fastSetPerThread(numThreads, fast_set(0));
    vector<vector<int>> tempIntDoubleSizePerThread(numThreads);
    vector<int> candidates(m_AdjacencyArray.size());
    vector<int> toRemove(m_AdjacencyArray.size());
    vector<char> isCandidate = vector<char>(m_AdjacencyArray.size());
    vector<unique_ptr<vector<int>>> tempBuffer(numThreads);
    vector<vector<bool>> markedVerticesPerThread(numThreads);
    for(int threadId = 0; threadId < numThreads; threadId++) {
        vector<int> tempInt1 = vector<int>(m_AdjacencyArray.size());
        tempInt1PerThread[threadId] = tempInt1;
        vector<int> tempInt2 = vector<int>(m_AdjacencyArray.size());
        tempInt2PerThread[threadId] = tempInt2;
        fast_set fastSet(m_AdjacencyArray.size());
        fastSetPerThread[threadId] = fastSet;
        vector<int> tempIntDoubleSize = vector<int>(m_AdjacencyArray.size() * 2);
        tempIntDoubleSizePerThread[threadId] = tempIntDoubleSize;
        tempBuffer[threadId] = unique_ptr<vector<int>>(new vector<int>());
        tempBuffer[threadId]->reserve(maxTempSize);
        markedVerticesPerThread[threadId] = vector<bool>(m_AdjacencyArray.size(), false);
    }

    int numUnconfinedReductions = 0;
    int numLPReductions = 0;
    
    assert(checkDegrees());

    double startClock = omp_get_wtime();
    double tmpClock;
    double LPTime = 0;
    double unconfinedTime = 0;
    double isolatedCliqueTime = 0;

    int numIterations = 0;
    while(true) {
        numIterations++;

        tmpClock = omp_get_wtime();
        bool changedIsolatedClique = RemoveIsolatedClique(markedVerticesPerThread, toRemove, tempBuffer);
        isolatedCliqueTime += omp_get_wtime() - tmpClock;
        if(changedIsolatedClique) continue;

        tmpClock = omp_get_wtime();
        bool changedUnconfined = RemoveUnconfined(fastSetPerThread, tempInt1PerThread, tempInt2PerThread, tempIntDoubleSizePerThread, isCandidate, candidates, toRemove, tempBuffer);
        unconfinedTime += omp_get_wtime() - tmpClock;
        if(changedUnconfined) continue;

        tmpClock = omp_get_wtime();
        bool changedLP = LPReduction();
        LPTime = omp_get_wtime() - tmpClock;
        if(changedLP) continue;
        break;
    }
    std::cout << "Num iterations: " << numIterations << std::endl;


    double endClock = omp_get_wtime();


    cout << "Number of vertices removed by unconfined reduction: " << numUnconfinedReductions << endl;
    cout << "Number of vertices removed by LP reduction: " << numLPReductions << endl;
    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;

    std::cout << "LP time: " << LPTime << std::endl;
    std::cout << "Unconfined time: " << unconfinedTime << std::endl;

    assert(checkDegrees());
}

bool parallel_reductions_fine_grained::checkDegrees() {
    for(int vertex = 0; vertex < neighbors.size(); ++vertex) if(inGraph.Contains(vertex)) {
        int deg = 0;
        for(int const neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            deg++;
        }
        if(vertexDegree[vertex] != deg) {
            std::cout << "Vertex " << vertex << " has degree " << deg << " but vertexDegree[vertex] has value " << vertexDegree[vertex] << " number of elements in neighbors[vertex]: " << neighbors[vertex].Size() << std::endl;
            return false;
        }
    }
    return true;
}

void parallel_reductions_fine_grained::ApplyKernelSolutionToReductions(vector<Reduction> const &vReductions)
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
