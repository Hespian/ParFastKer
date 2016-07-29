#include <stdio.h>
#include <string.h>
#include <iostream>
#include <argtable2.h>
#include <dirent.h>

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
#include <iostream>
#include <fstream>

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

    auto dir = opendir(partitions_directory.c_str());
    std::vector<std::string> partition_files;
    struct dirent *dp;
    while ((dp = readdir(dir)) != NULL)
       if(ends_with(dp->d_name, ".partition")) {
        partition_files.push_back(dp->d_name);
       }
    (void)closedir(dir);

    std::vector<unsigned int> vertexTimes = std::vector<unsigned int>(G.number_of_nodes(),0);
    std::cout << "INIT!!!" << std::endl;

    for(std::string partition_file: partition_files) {
        std::cout << "---------------------------------------------------------------------" << std::endl;
        std::cout << "Number of blocks: " << partition_file.substr(0, partition_file.find ('.')) << std::endl;
        std::string partition_file_path = "";
        partition_file_path += partitions_directory;
        partition_file_path += "/";
        partition_file_path += partition_file;
        graph_io::readPartition(G, partition_file_path);
        std::vector<int> partitions(G.number_of_nodes());
        forall_nodes(G, node) {
            partitions[node] = G.getPartitionIndex(node);
        } endfor

        for(int i = 0; i < mis_config.num_reps; ++i) {
            std::cout << "---------------------------------------------------------------------" << std::endl;
            std::cout << "New repitition: " << i  << std::endl;
            std::unique_ptr<full_reductions> full_reducer_parallel = std::unique_ptr<full_reductions>(new full_reductions(adj_for_parallel_aglorithm, partitions));

            full_reducer_parallel->reduce_graph(vertexTimes);
        }
    }

    std::ofstream outputFile(mis_config.graph_filename + "-workload_weights.weights");
    if (outputFile.is_open()) {
        for(int vertex = 0; vertex < G.number_of_nodes(); vertex++) {
            outputFile << vertexTimes[vertex] + 1 << "\n";
        }
        outputFile.close();
    }
    else { 
        std::cout << "Unable to open file";
        exit(1);
    }


    return 0;
}
