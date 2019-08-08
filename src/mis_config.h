/**
 * mis_config.h
 * Purpose: Configuration used for the evolutionary maximum independent set algorithms.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
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

#ifndef _MIS_CONFIG_H_
#define _MIS_CONFIG_H_

#include <string>

#include "definitions.h"

// Configuration for the calculation of the MIS
struct MISConfig {
    // Name of the graph file.
    std::string graph_filename;
    // Directory containing partitions
    std::string partition_directory;
    // Name of the output file.
    std::string output_filename;
    // Number of repititions for benchmark
    int num_reps;
    // Write the log into a file
    bool print_log;
    // Write the inpendent set into a file
    bool write_graph;
    // Write the log into the console
    bool console_log;
    // Apply all reductions to reduce the graph size
    bool all_reductions;
    // Check graph sortedness
    bool check_sorted;
};

#endif
