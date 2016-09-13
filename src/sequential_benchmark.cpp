#include <stdio.h>
#include <string.h>
#include <iostream>
#include <argtable2.h>
#include <dirent.h>
#include <omp.h>

#include "timer.h"
// #include "ils/ils.h"
// #include "ils/local_search.h"
#include "mis_log.h"
#include "graph_io.h"
// #include "reduction_evolution.h"
#include "mis_config.h"
// #include "greedy_mis.h"
#include "parse_parameters.h"
#include "data_structure/graph_access.h"
// #include "data_structure/mis_permutation.h"
#include "sequential/branch_and_reduce_algorithm.h"
#include "full_reductions.h"
#include <memory>

inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int main(int argn, char **argv) {
    mis_log::instance()->print_title();
    
    MISConfig mis_config;
    std::string graph_filepath;
    std::string partitions_directory;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath, partitions_directory);
    if (ret_code) {
        return 0;
    }
    mis_config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of( '/' ) +1);
    mis_log::instance()->set_config(mis_config);

    // Read the graph
    graph_access G;
    graph_io::readGraphWeighted(G, graph_filepath);
    mis_log::instance()->set_graph(G);
    
    // Print setup information
    mis_log::instance()->print_graph();
    mis_log::instance()->print_config();

    // initialize full reducer
    std::cout << "Creating graph" << std::endl;
    std::vector<std::vector<int>> adj_for_parallel_aglorithm(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
        adj_for_parallel_aglorithm[node].reserve(G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
            adj_for_parallel_aglorithm[node].push_back(neighbor);
        } endfor
    } endfor
    std::cout << "Finished creating graph" << std::endl;

    std::vector<int> partitions(G.number_of_nodes());
    forall_nodes(G, node) {
        partitions[node] = 0;
    } endfor

    omp_set_num_threads(1);
    for(int i = 0; i < mis_config.num_reps; ++i) {
        std::cout << "---------------------------------------------------------------------" << std::endl;
        std::cout << "New repitition: " << i  << std::endl;
        std::unique_ptr<parallel_reductions> reducer = std::unique_ptr<parallel_reductions>(new parallel_reductions(adj_for_parallel_aglorithm, partitions));
        reducer->reduce_graph_parallel();
        std::cout << "Kernel size: " << reducer->size() << std::endl;
    }


    return 0;
}
