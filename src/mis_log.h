/******************************************************************************
 * mis_log.h
 *
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

#ifndef _MIS_LOG_H_
#define _MIS_LOG_H_

#include <sstream>

#include "timer.h"
#include "mis_config.h"
#include "data_structure/graph_access.h"

class mis_log {
    public:
        /**
         * Get the singleton logger instance.
         * 
         * @return Instance of the logger.
         */
        static mis_log *instance() {
            static mis_log inst;
            return &inst;
        };

        /**
         * Set the config.
         *
         * @param config Config for the evolutionary algorithm.
         */
        void set_config(MISConfig & config);

        /**
         * Set the graph.
         *
         * @param G Graph representation.
         */
        void set_graph(graph_access & G);

        /**
         * Add a newline to the log.
         */
        void print_newline();

        /**
         * Print information about the graph.
         */
        void print_graph();

        /**
         * Print the current config.
         */
        void print_config();

        /**
         * Print information about a round.
         *
         * @param mis_config Config for the logger.
         */
        // void print_round(MISConfig & mis_config);

        /**
         * Print information about the reduction step.
         * Includes number of extracted nodes and resulting kernel size
         *
         * @param mis_config Config for the logger.
         * @param extracted_nodes Number of removed nodes.
         * @param kernel_size Number of remaining nodes.
         */
        void print_reduction(MISConfig & mis_config, unsigned int extracted_nodes, unsigned int kernel_size);

        /**
         * Print information about a single repetition.
         *
         * @param mis_config Config for the logger.
         */
        // void print_repetition(MISConfig & mis_config);

        /** 
         * Print the time needed to construct the separator pool.
         */
        // void print_separator();

        /**
         * Print banner for the partitioning phase.
         */
        // void print_pool_title();

        /**
         * Print banner for the evolution phase.
         */
        // void print_evolution_title();

        /**
         * Print banner for the initialization phase.
         */
        // void print_init_title();

        /**
         * Print the final results.
         */
        // void print_results();

        /**
         * Print a title.
         */
        void print_title();
        
        /**
         * Restart the timer for the total time including IO, etc.
         */
        // void restart_total_timer();

        /**
         * Restart the timer for the evolutionary algorithm.
         */
        // void restart_evo_timer();

        /**
         * Get the timer for the evolutionary algorithm.
         */
        // double get_evo_timer();

        /**
         * Restart the timer for a single operation.
         */
        // void restart_operator_timer();

        /**
         * Restart the timer for the separator pool.
         */
        // void restart_building_pool_timer();

        /**
         * Time taken for building the separator pool.
         * Restarts the after pool timer.
         *
         * @return Time taken.
         */
        // double get_pool_building_time();

        /**
         * Time taken after building the separator pool.
         *
         * @return Time taken.
         */
        // double get_after_pool_time();
    
        /**
         * Increment the number of rounds by one.
         */
        // void inc_rounds();

        /**
         * Get current rount.
         * 
         * @return Current round of the evolutionary algorithm.
         */
        // unsigned int current_round() { return number_of_rounds; };

        /** 
         * Increment the number of repetitions by one.
         */
        // void inc_repetitions();

        /**
         * Get current repetition.
         *
         * @return Current repetition of the evolutionary algorithm.
         */
        // unsigned int current_repetition() { return number_of_repetitions; };

        /**
         * Set the current operator.
         *
         * @param operator_name Name of the operator.
         */
        // void set_operator(std::string operator_name);

        /**
         * Set the size of the result of the operator.
         *
         * @param result Solution size.
         */
        // void set_result_operator(unsigned int result);

        /**
         * Update the size of the best solution.
         *
         * @param mis_config Config for the logger.
         * @param size Candidate to replace the best solution size.
         */
        // void set_best_size(MISConfig & mis_config, unsigned int size);

        /**
         * Reset the size of the best solution.
         */
        // void reset_best_size();

        /**
         * Set the average solution size.
         *
         * @param avg_size Average solution size.
         */
        // void set_avg_solution_size(double avg_size);


    private:
        // General information
        timer total_timer;
        timer evo_timer;
        timer operator_timer;
        timer pool_timer;
        timer non_pool_timer;
        std::stringstream filebuffer_string;
        MISConfig log_config;

        // Graph informations
        std::string graph_name;
        unsigned int number_of_nodes;
        unsigned int number_of_edges;
        unsigned int arc_scans;
        double avg_degree;
        double density;
        
        // Evolutionary informations
        unsigned int number_of_rounds;
        unsigned int number_of_repetitions;
        unsigned int best_solution_size;
        double avg_solution_size;

        // Repetition information
        std::string evo_operator;
        unsigned int result_operator;

        // Reduction information
        unsigned int current_is_size;
        unsigned int optimum_size;

        // Results information
        double total_time_taken;

        unsigned int separator_selected;
        unsigned int cover_selected;
        unsigned int multiway_selected;

        unsigned int separator_improv;
        unsigned int cover_improv;
        unsigned int multiway_improv;

        double total_separator_time;
        double total_cover_time;
        double total_multiway_time;

        double avg_separator_time;
        double avg_cover_time;
        double avg_multiway_time;

        unsigned int total_separator_improv;
        unsigned int total_cover_improv;
        unsigned int total_multiway_improv;

        double avg_separator_improv;
        double avg_cover_improv;
        double avg_multiway_improv;

        double time_taken_best;
        std::string operator_best;
        unsigned int repetition_best;
        unsigned int round_best;

        double time_for_building_pool;
        double time_since_building_pool;

        /**
         * Default Constructor.
         */
        mis_log();

        /**
         * Default Destructor.
         */
        virtual ~mis_log();

        /**
         * Compute the average time for each operation.
         */
        void compute_avg();
};

#endif
