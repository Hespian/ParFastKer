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
#include <limits>

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

    auto is_base = full_reducer_parallel->get_current_is_size_with_folds();
    mis_log::instance()->print_reduction(mis_config, is_base, full_reducer_parallel->number_of_nodes_remaining());


    std::cout << "Generating new graph for sequential solver" << std::endl;
    std::vector<int> parallel_to_sequential_map(G.number_of_nodes());
    int nodecount = 0;
    forall_nodes(G, node) {
    	if(full_reducer_parallel->x[node] < 0) {
	        parallel_to_sequential_map[node] = nodecount++;
	    }
    } endfor



    // initialize full reducer
    std::vector<std::vector<int>> adj_for_sequential_aglorithm(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
    	if(full_reducer_parallel->x[node] < 0) {
	        adj_for_sequential_aglorithm[parallel_to_sequential_map[node]].reserve(G.getNodeDegree(node));
	        for(auto neighbor : full_reducer_parallel->adj[node]) {
	            if(full_reducer_parallel->x[neighbor] < 0) {
	            	adj_for_sequential_aglorithm[parallel_to_sequential_map[node]].push_back(parallel_to_sequential_map[neighbor]);
	            }
	        }
	    }
    } endfor

    std::cout << "Solving kernel with sequential algorithm" << std::endl;

    std::unique_ptr<branch_and_reduce_algorithm> full_reducer_sequential = std::unique_ptr<branch_and_reduce_algorithm>(new branch_and_reduce_algorithm(adj_for_sequential_aglorithm, adj_for_sequential_aglorithm.size()));

    timer t = timer();
    double timelimit = std::numeric_limits<double>::max();
    full_reducer_sequential->solve(t, timelimit);



    std::cout << "Applying solution to parallel result" << std::endl;
    forall_nodes(G, node) {
    	if(full_reducer_parallel->x[node] < 0) {
	        full_reducer_parallel->x[node] = full_reducer_sequential->y[parallel_to_sequential_map[node]];
	    }
    } endfor

    full_reducer_parallel->undoReductions();

    std::cout <<  "checking solution validity ..."  << std::endl;

    int counter = 0;
    forall_nodes(G, node) {
            if( full_reducer_parallel->x[node] == 0 ) {
                    counter++;
                    forall_out_edges(G, e, node) {
                            NodeID target = G.getEdgeTarget(e);
                            if(full_reducer_parallel->x[target] == 0) {
                                std::cout <<  "not an independent set! "  << node << " and " << target << " are both in the is!" << std::endl;
                                exit(1);
                            }
                    } endfor
            }
    } endfor
    std::cout <<  "valid"  << std::endl; 

    std::cout <<  "checking solution size ..."  << std::endl;
    std::cout <<  "Solving without parallel reductions ..."  << std::endl;
    // initialize full reducer
    std::vector<std::vector<int>> adj_for_sequential_aglorithm_test(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
        adj_for_sequential_aglorithm_test[node].reserve(G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
        	adj_for_sequential_aglorithm_test[node].push_back(neighbor);
        } endfor
    } endfor

    std::unique_ptr<branch_and_reduce_algorithm> full_reducer_sequential_test = std::unique_ptr<branch_and_reduce_algorithm>(new branch_and_reduce_algorithm(adj_for_sequential_aglorithm_test, adj_for_sequential_aglorithm_test.size()));

    t = timer();
    full_reducer_sequential_test->solve(t, timelimit);

    std::cout <<  "checking sequential solution validity ..."  << std::endl;
    forall_nodes(G, node) {
            if( full_reducer_sequential_test->y[node] == 0 ) {
                    counter++;
                    forall_out_edges(G, e, node) {
                            NodeID target = G.getEdgeTarget(e);
                            if(full_reducer_sequential_test->y[target] == 0) {
                                std::cout <<  "Not an independent set! "  << node << " and " << target << " are both in the is!" << std::endl;
                                exit(1);
                            }
                    } endfor
            }
    } endfor
    std::cout <<  "valid"  << std::endl; 

    int sequential_size = 0;
    for(auto i : full_reducer_sequential_test->y) {
    	if(i == 0) {
    		sequential_size++;
    	}
    }

    int parallel_size = 0;
    for(auto i : full_reducer_parallel->x) {
    	if(i == 0) {
    		parallel_size++;
    	}
    }

    if(sequential_size != parallel_size){
        std::cout << "Sequential size: " << sequential_size << "; With parallel kernel: " << parallel_size << std::endl;
        exit(1);
    }

    std::cout << "Same size" << std::endl;

    return 0;
}