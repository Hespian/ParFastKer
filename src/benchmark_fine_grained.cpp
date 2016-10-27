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
#include "parallel_reductions_fine_grained.h"
#include <memory>

inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool testCorrectness(std::unique_ptr<parallel_reductions_fine_grained> &full_reducer_parallel, graph_access &G) {

    std::vector<std::vector<int>> adj_for_sequential_aglorithm = full_reducer_parallel->getKernel();

    std::cout << "Solving kernel with sequential algorithm" << std::endl;

    std::unique_ptr<branch_and_reduce_algorithm> full_reducer_sequential = std::unique_ptr<branch_and_reduce_algorithm>(new branch_and_reduce_algorithm(adj_for_sequential_aglorithm, adj_for_sequential_aglorithm.size()));

    timer t = timer();
    double timelimit = std::numeric_limits<double>::max();
    full_reducer_sequential->solve(t, timelimit);



    std::cout << "Applying solution to parallel result" << std::endl;
    full_reducer_parallel->applyKernelSolution(full_reducer_sequential->y);

    std::cout <<  "checking solution validity ..."  << std::endl;
    std::vector<int> parallel_solution = full_reducer_parallel->independent_set;

    int counter = 0;
    forall_nodes(G, node) {
            if( parallel_solution[node] == 0 ) {
                    counter++;
                    forall_out_edges(G, e, node) {
                            NodeID target = G.getEdgeTarget(e);
                            if(parallel_solution[target] == 0) {
                                std::cout <<  "not an independent set! "  << node << " and " << target << " are both in the is!" << std::endl;
                                exit(1);
                            }
                    } endfor
            } else if (parallel_solution[node] < 0) {
                std::cout << "Incomplete Solution! Vertex " << node << " undefined!" << std::endl;
                exit(1);
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
            } else if (full_reducer_sequential_test->y[node] == -1) {
                std::cout <<  "Incomplete solution! "  << node <<  " not determined!" << std::endl;
                exit(1);
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
    for(auto i : parallel_solution) {
        if(i == 0) {
            parallel_size++;
        }
    }

    if(sequential_size != parallel_size){
        std::cout << "Sequential size: " << sequential_size << "; With parallel kernel: " << parallel_size << std::endl;
        exit(1);
    }

    std::cout << "Same size" << std::endl;

    return true;
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

    omp_set_num_threads(4);
    for(int i = 0; i < mis_config.num_reps; ++i) {
        std::cout << "---------------------------------------------------------------------" << std::endl;
        std::cout << "New repitition: " << i  << std::endl;
        std::unique_ptr<parallel_reductions_fine_grained> reducer = std::unique_ptr<parallel_reductions_fine_grained>(new parallel_reductions_fine_grained(adj_for_parallel_aglorithm));
        reducer->reduce_graph_parallel();
        std::cout << "Kernel size: " << reducer->size() << std::endl;
        // testCorrectness(reducer, G);
    }


    return 0;
}
