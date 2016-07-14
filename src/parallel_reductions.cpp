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
#define REMOVE_NEIGHBOR(partition, neighbor, vertex) if(partition == partitions[neighbor]) neighbors[neighbor].Remove(vertex); else neighborhoodChanged.Insert(neighbor);

using namespace std;

ProfilingHelper_t profilingHelper;

parallel_reductions::parallel_reductions(vector<vector<int>> const &adjacencyArray)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size(), true)
 , neighborhoodChanged(adjacencyArray.size(), false)
 , partitions(adjacencyArray.size())
#ifdef TIMERS
 , replaceTimer(0)
 #endif // TIMERS
 , m_bAllowVertexFolds(true)
{
    std::cout << "Start constructor" << std::endl;
    for (size_t u=0; u < adjacencyArray.size(); ++u) {
        neighbors[u].InitializeFromAdjacencyArray(m_AdjacencyArray, u);
    }
    independent_set.resize(m_AdjacencyArray.size());
    std::cout << "Finished constructor" << std::endl;
}

parallel_reductions::~parallel_reductions()
{

#ifdef TIMERS
    cout << "Total time spent undoing  reductions  : " << (replaceTimer/(double)CLOCKS_PER_SEC) << endl;
#endif // TIMERS
}

void parallel_reductions::partitionGraph(int numPartitions, string partitioner) {
    std::cout << "Start partitioning into " << numPartitions << " partitions using " << partitioner << std::endl;
    partition_nodes.resize(numPartitions);
    int N = m_AdjacencyArray.size();
    if (partitioner == "kahip") {
        std::vector<int> xadj;
        std::vector<int> adjncy;
        for(auto node_adj_list : m_AdjacencyArray) {
            xadj.push_back(adjncy.size());
            for(auto neighbor : node_adj_list) {
                adjncy.push_back(neighbor);
            }
        }
        xadj.push_back(adjncy.size());
        int edgecut = 0;
        int number_of_partitions = numPartitions;
        double imbalance = 0.03;
        kaffpa(&N, NULL, xadj.data(), NULL, adjncy.data(), &number_of_partitions, &imbalance, false, 1337, FASTSOCIAL, &edgecut, partitions.data());
        std::cout << "Edgecut: " << edgecut << std::endl;
        for(int node = 0; node < N; ++node) {
            partition_nodes[partitions[node]].push_back(node);
        }
    } else {


        int edgecount = 0;
        for(auto node_adj : m_AdjacencyArray) {
            edgecount += node_adj.size();
        }
        edgecount /= 2;

        ofstream outputFile("tmpgraph.graph");
        if (outputFile.is_open()) {
            outputFile << N << " " << edgecount << " 0\n";
            for(auto node_adj : m_AdjacencyArray) {
                for(auto neighbor : node_adj) {
                    outputFile << neighbor + 1 << " ";
                }
                outputFile << "\n";
            }
            outputFile.close();
        }
        else { 
            cout << "Unable to open file";
            exit(1);
        }

        std::ostringstream oss;
        if (partitioner == "parallel_kahip")
            oss << "mpirun -n " << numPartitions << " ../../parallel_social_partitioning_package/deploy/parallel_label_compress ./tmpgraph.graph --k=" << numPartitions << " --preconfiguration=ultrafast --seed 1337 >> partition_output";
        else if (partitioner == "lpa")
            oss << "../../KaHIPLPkway/deploy/label_propagation --k " << numPartitions << " ./tmpgraph.graph --seed=1337 --label_propagation_iterations=1 --output_filename=tmppartition >> partition_output";
        else {
            cout << "Unknown partitioner" << std::endl;
            exit(1);
        }
        std::string command = oss.str();
        std::cout << command << std::endl;
        int system_succesfull = system(command.c_str());
        system("rm tmpgraph.graph");
        if(system_succesfull != 0) {
            cout << "Command unsuccessful" << std::endl;
            exit(1);
        }

        std::string line;
        std::ifstream inputfile("./tmppartition");
        if (!inputfile) {
                std::cerr << "Error opening file" << std::endl;
                exit(1);
        }
        for(int node  = 0; node < N; ++node) {
            std::getline(inputfile, line);
                if (line[0] == '%') { //Comment
                        node--;
                        continue;
                }
                partitions[node] = atol(line.c_str());
        }

        inputfile.close();
        system("rm tmppartition");

        for(int node = 0; node < N; ++node) {
            partition_nodes[partitions[node]].push_back(node);
        }
    }
    std::cout << "Finished partitioning" << std::endl;
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
    for(auto Reductions: ReductionsPerPartition) {
        ApplyKernelSolutionToReductions(Reductions);
    }
}

bool parallel_reductions::RemoveIsolatedClique(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &isolatedCliqueCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    /*if(neighbors[vertex].Size() > ISOLATED_CLIQUE_MAX_NEIGHBORS)
        return false;*/

    for (int const neighbor : neighbors[vertex]) {
        if (partitions[neighbor] != partition) {
            profilingAddTimeUnsuccessfulIsolatedCliqueDegreeOrPartition(&profilingHelper, partition);
            return false;
        }
        if (neighbors[neighbor].Size() < neighbors[vertex].Size()) {
            profilingAddTimeUnsuccessfulIsolatedCliqueDegreeOrPartition(&profilingHelper, partition);
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
            remaining.Remove(neighbor);
            independent_set[neighbor] = 1;
            // reduction.AddNeighbor(neighbor);
            // reduction.AddRemovedEdge(vertex,   neighbor);
            // reduction.AddRemovedEdge(neighbor, vertex);
        }
        inGraph.Remove(vertex);

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

bool parallel_reductions::isTwoNeighborhoodInSamePartition(int const vertex, int const partition) {
    for(int neighbor1 : neighbors[vertex]) {
        if(partitions[neighbor1] != partition) {
            return false;       
        }
        updateNeighborhood(neighbor1);
        for(int neighbor2 : neighbors[neighbor1]) {
            if(partitions[neighbor2] != partition) {
                return false;       
            }
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
    if (neighbors[neighbors[vertex][0]].Contains(neighbors[vertex][1])) {
        profilingAddTimeUnsuccessfulFoldAdjacent(&profilingHelper, partition);
        return false; // neighbors can't be adjacent.
    }
    if (!isTwoNeighborhoodInSamePartition(vertex, partition)) { 
        profilingAddTimeUnsuccessfulFoldWrongPartition(&profilingHelper, partition);
        return false;
    }

    foldedVertexCount++;

    int const vertex1(neighbors[vertex][0]);
    int const vertex2(neighbors[vertex][1]);

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
    for(int neighbor: neighbors[vertex]) {
        if(!inGraph.Contains(neighbor)) {
            verticesToRemove.push_back(neighbor);
        }
    }
    for(int neighborToRemove : verticesToRemove) {
        neighbors[vertex].Remove(neighborToRemove);
    }

    profilingAddTimeUpdateNeighborhood(&profilingHelper, partitions[vertex]);
}

void parallel_reductions::reduce_graph(int numPartitions, string partitioner) {
    profilingInit(&profilingHelper, &neighbors, numPartitions);
    partitionGraph(numPartitions, partitioner);

    vector<vector<bool>> vMarkedVerticesPerPartition(numPartitions);
    vector<ArraySet> remainingPerPartition(numPartitions);
    ReductionsPerPartition = vector<vector<Reduction>>(numPartitions);
    for(int partition = 0; partition < numPartitions; partition++) {
        std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
        vMarkedVerticesPerPartition[partition] = vMarkedVertices;
        ArraySet remaining = ArraySet(m_AdjacencyArray.size());
        remainingPerPartition[partition] = remaining;
    }
    
    double startClock = omp_get_wtime();

    #pragma omp parallel for
    for(int partition = 0; partition < numPartitions; partition++) {
        ApplyReductions(partition, partition_nodes[partition], ReductionsPerPartition[partition], vMarkedVerticesPerPartition[partition], remainingPerPartition[partition]);
    }


    double endClock = omp_get_wtime();
    profilingPrint(&profilingHelper);
    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
}

void parallel_reductions::ApplyReductions(int const partition, vector<int> vertices, vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining)
{
    double startClock = omp_get_wtime();
    std::cout << partition << ": Filling remaining vertices set..." << std::endl;
    remaining.Clear();
    for (int const vertex : vertices) {
        INSERT_REMAINING(partition, remaining, vertex);
    }
    int iterations(0);
    int foldedVertexCount(0);
    int isolatedCliqueCount(0);
    std::cout << partition << ": Starting reductions..." << std::endl;
    while (!remaining.Empty()) {
        int const vertex = *(remaining.begin());
        remaining.Remove(vertex);
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
    cout << partition << ": Time spent applying reductions  : " << (endClock - startClock) << endl;
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
