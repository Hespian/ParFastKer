#include <stdio.h>
#include <string.h>
#include <iostream>
#include <argtable2.h>

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
#include "parallel_branch_and_reduce_algorithm.h"
#include <memory>

int main(int argn, char **argv) {
    mis_log::instance()->restart_total_timer();
    mis_log::instance()->print_title();
    
    MISConfig mis_config;
    std::string graph_filepath;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
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
    std::vector<std::vector<int>> adj_for_sequential_aglorithm(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
        adj_for_sequential_aglorithm[node].reserve(G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
            adj_for_sequential_aglorithm[node].push_back(neighbor);
        } endfor
    } endfor

    std::unique_ptr<branch_and_reduce_algorithm> full_reducer_sequential = std::unique_ptr<branch_and_reduce_algorithm>(new branch_and_reduce_algorithm(adj_for_sequential_aglorithm, adj_for_sequential_aglorithm.size()));

    full_reducer_sequential->reduce_graph();

    size_t is_base = full_reducer_sequential->get_current_is_size_with_folds();
    mis_log::instance()->print_reduction(mis_config, is_base, full_reducer_sequential->number_of_nodes_remaining());    

    // initialize full reducer
    std::vector<std::vector<int>> adj_for_parallel_aglorithm(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
        adj_for_parallel_aglorithm[node].reserve(G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
            adj_for_parallel_aglorithm[node].push_back(neighbor);
        } endfor
    } endfor

    std::unique_ptr<parallel_branch_and_reduce_algorithm> full_reducer_parallel = std::unique_ptr<parallel_branch_and_reduce_algorithm>(new parallel_branch_and_reduce_algorithm(adj_for_parallel_aglorithm, adj_for_parallel_aglorithm.size(), mis_config));

    full_reducer_parallel->reduce_graph();

    is_base = full_reducer_parallel->get_current_is_size_with_folds();
    mis_log::instance()->print_reduction(mis_config, is_base, full_reducer_parallel->number_of_nodes_remaining());



   /* std::cout <<  "checking solution ..."  << std::endl;
    // initialize full reducer
    std::vector<std::vector<int>> adj_for_parallel_aglorithm_test(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
        adj_for_parallel_aglorithm_test[node].reserve(G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
            adj_for_parallel_aglorithm_test[node].push_back(neighbor);
        } endfor
    } endfor

    std::unique_ptr<parallel_branch_and_reduce_algorithm> full_reducer_parallel_for_test = std::unique_ptr<parallel_branch_and_reduce_algorithm>(new parallel_branch_and_reduce_algorithm(adj_for_parallel_aglorithm_test, adj_for_parallel_aglorithm_test.size(), mis_config));

    std::vector<int> independent_set = full_reducer_parallel_for_test->compute_maximal_is();

    int counter = 0;
    forall_nodes(G, node) {
            if( independent_set[node] == 0 ) {
                    counter++;
                    forall_out_edges(G, e, node) {
                            NodeID target = G.getEdgeTarget(e);
                            if(independent_set[target] == 0) {
                                std::cout <<  "not an independent set!"  << std::endl;
                                exit(0);
                            }
                    } endfor
            }
    } endfor
    std::cout <<  "done ..."  << std::endl; */
    return 0;
}