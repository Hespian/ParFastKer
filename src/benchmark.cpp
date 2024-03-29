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
#include "omp.h"
#include <algorithm>

inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline bool starts_with(std::string const & value, std::string const & start)
{
  if (start.size() > value.size()) return false;
  return std::equal(start.begin(), start.end(), value.begin());
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

    // auto outerDir=opendir(partitions_directory.c_str());
    // std::vector<std::string> weight_dirs;
    // struct dirent *dp1;
    // while((dp1 = readdir(outerDir)) != NULL) { 
    //   if(starts_with(dp1->d_name, "weight")) {
	// weight_dirs.push_back(dp1->d_name);
    //   }
    // }
    // (void)closedir(outerDir);

    std::string partitions_dir_original = std::string(partitions_directory);
    
    // for(std::string weight_dir : weight_dirs) {
    //   if(!ends_with(weight_dir, "one_ultrafast")) continue;
    //   //if(!ends_with(weight_dir, "LPA")) continue;
    //   std::cout << "========================================================================" << std::endl;
    //   std::cout << weight_dir << std::endl;
    //   partitions_directory = partitions_dir_original + "/" + weight_dir;


    // auto dir = opendir(partitions_directory.c_str());
    std::vector<std::string> partition_files;
    // struct dirent *dp;
    // while ((dp = readdir(dir)) != NULL)
    //    if(ends_with(dp->d_name, ".partition")) {
    //     partition_files.push_back(dp->d_name);
    //    }
    // (void)closedir(dir);

    // int max_blocks = 0;
    // for(std::string partition_file:partition_files) {
    //   int numPartitions = std::stoi(partition_file.substr(0, partition_file.find ('.')));
    //   if(numPartitions > max_blocks) {
	// max_blocks = numPartitions;
    //   }
    // }

    partition_files.push_back(partitions_dir_original);

    for(std::string partition_file: partition_files) {
        // int numPartitions = std::stoi(partition_file.substr(12, partition_file.length()));
        std::cout << partition_file.substr(partition_file.rfind ('/') + 1, partition_file.rfind ('.') - partition_file.rfind ('/')) << std::endl;
        int numPartitions = std::stoi(partition_file.substr(partition_file.rfind ('/') + 1, partition_file.rfind ('.') - partition_file.rfind ('/')));
        int max_threads = 32;
        omp_set_num_threads(std::min(numPartitions, max_threads));
        if(numPartitions == 1) {
            omp_set_num_threads(2);
        }
        if(numPartitions > 32)
            continue;
        // if(numPartitions != 32) {
        //     continue;
        // }
        // if(numPartitions != 1) {
        //     continue;
        // }
        //if(numPartitions != max_blocks)
        //  continue;
        std::cout << "---------------------------------------------------------------------" << std::endl;
        std::cout << "Number of blocks: " << numPartitions << std::endl;
        std::string partition_file_path = "";
        partition_file_path += partitions_directory;
        // partition_file_path += "/";
        // partition_file_path += partition_file;
        graph_io::readPartition(G, partition_file_path);
        std::vector<int> partitions(G.number_of_nodes());
        forall_nodes(G, node) {
            partitions[node] = G.getPartitionIndex(node);
        } endfor

              for(int i = 0; i < mis_config.num_reps; ++i) {
                  std::cout << "---------------------------------------------------------------------" << std::endl;
                  std::cout << "New repitition: " << i  << std::endl;
                  std::unique_ptr<full_reductions> full_reducer_parallel = std::unique_ptr<full_reductions>(new full_reductions(adj_for_parallel_aglorithm, partitions));

                  bool writeKernel = false; //(i == mis_config.num_reps - 1 && numPartitions == 32);
                  // std::string kernel_path = "kernels/" + mis_config.graph_filename + ".quasikernel";
                  std::string kernel_path = mis_config.output_filename;
                  full_reducer_parallel->reduce_graph(mis_config.write_graph, kernel_path);

                  // if(i == mis_config.num_reps - 1 && numPartitions == 32) {
                  //     int numKernelEdges = 0;
                  //     std::vector<std::vector<int>> kernel = full_reducer_parallel->getKernel();
                  //     for(std::vector<int> neighbors : kernel) {
                  //         numKernelEdges += neighbors.size();
                  //     }
                  //     std::string kernel_path = "kernels/" + mis_config.graph_filename + ".kernel";
                  //     std::cout << "Writing kernel to " << kernel_path << std::endl;
                  //     std::ofstream f(kernel_path.c_str());
                  //     f << kernel.size() <<  " " <<  numKernelEdges / 2 << std::endl;

                  //     for(auto vertexNeighbors : kernel) {
                  //         std::sort(vertexNeighbors.begin(), vertexNeighbors.end());
                  //         for(auto neighbor : vertexNeighbors) {
                  //             f <<   neighbor + 1 << " " ;
                  //         } 
                  //         f <<  std::endl;
                  //     }

                  //     f.close();
                  // }
              }
    }
    // }


    return 0;
}
