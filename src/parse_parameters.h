/**
 * parse_parameters.h
 * Purpose: Parse command line parameters.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <omp.h>

#include "configuration_mis.h"

/**
 * Parse the given parameters and apply them to the config.
 *
 * @param argn Number of parameters.
 * @param argv Values of the parameters.
 * @param mis_config Config to store the values in.
 * @param graph_filename String to store the filename of the graph.
 *
 * @return -1 if there was an error. 0 otherwise.
 */
int parse_parameters(int argn, char **argv,
                     MISConfig & mis_config,
                     std::string & graph_filename) {
    const char *progname = argv[0];

    // Setup the argtable structs
    struct arg_lit *help                = arg_lit0(NULL, "help", "Print help.");
    struct arg_int *user_seed           = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
    struct arg_str *user_conf           = arg_str0(NULL, "config", NULL, "Configuration to use. ([standard|social|full_standard|full_social]). Standard/social use different modes of the graph partitioning tool. Full configurations use more time consuming parameters.");
    
    struct arg_str *partitioner         = arg_str0(NULL, "partitioner", NULL, "Partitioner to use. ([kahip, parallel_kahip, lpa]).");
    struct arg_int *kahip_mode          = arg_int0(NULL, "kahip_mode", NULL, "KaHIP mode to use.");
    struct arg_int *num_partitions      = arg_int0(NULL, "num_partitions", NULL, "Number of partitions to use.");

    struct arg_str *filename            = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
    struct arg_str *output              = arg_str0(NULL, "output", NULL, "Path to store resulting independent set.");
    struct arg_lit *console_log         = arg_lit0(NULL, "console_log", "Stream the log into the console");
    struct arg_lit *disable_checks      = arg_lit0(NULL, "disable_checks", "Disable sortedness check during I/O.");

    struct arg_end *end                 = arg_end(100);

    // Setup the argtable
    void *argtable[] = {
            help, 
            filename, 
            output,
            user_seed, 
            user_conf, 
            partitioner,
            kahip_mode,
            num_partitions,
            console_log,
            disable_checks,
            end
    };    

    // Choose standard configuration
    configuration_mis cfg;
    cfg.standard(mis_config);
    
    // Parse the arguments
    int nerrors = arg_parse(argn, argv, argtable);

    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable, "  %-40s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }

    if (nerrors > 0) {
        arg_print_errors(stderr, end, progname);
        printf("Try '%s --help' for more information.\n", progname);
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }

    if (user_conf->count > 0) {
        if (strcmp(user_conf->sval[0], "standard") == 0) cfg.standard(mis_config);
        else if (strcmp(user_conf->sval[0], "social") == 0) cfg.social(mis_config);
        else if (strcmp(user_conf->sval[0], "full_standard") == 0) cfg.full_standard(mis_config);
        else if (strcmp(user_conf->sval[0], "full_social") == 0) cfg.full_social(mis_config);
    }

    if (partitioner->count > 0) {
        if (strcmp(partitioner->sval[0], "kahip") == 0) mis_config.partitioner = "kahip";
        else if (strcmp(partitioner->sval[0], "parallel_kahip") == 0) mis_config.partitioner = "parallel_kahip";
        else if (strcmp(partitioner->sval[0], "lpa") == 0) mis_config.partitioner = "lpa";
    }

    if (filename->count > 0) {
        graph_filename = filename->sval[0];
    }   

    if (kahip_mode->count > 0) {
        mis_config.kahip_mode = kahip_mode->ival[0];
    }

    if (num_partitions->count > 0) {
        mis_config.number_of_partitions = num_partitions->ival[0];
    }

    if (user_seed->count > 0) {
        mis_config.seed = user_seed->ival[0];
    }

    // if (use_multiway_ns->count > 0) {
    //     mis_config.use_multiway_vc = false;
    // }

    // if (use_multiway_vc->count > 0) {
    //     mis_config.use_multiway_vc = true;
    // }

    // if (use_hopcroft->count > 0) {
    //     mis_config.use_hopcroft = true;
    // }

    // if (repetitions->count > 0) {
    //     mis_config.repetitions = repetitions->ival[0];
    // }

    //if (best_degree_frac->count > 0) {
        //mis_config.remove_fraction = best_degree_frac->dval[0];
    //}

    if (console_log->count > 0) {
        mis_config.console_log = true;
        mis_config.print_log = false;
    } else {
        mis_config.print_log = true;
    }

    if (disable_checks->count > 0) {
        mis_config.check_sorted = false;
    }

    if (output->count > 0) {
        mis_config.output_filename = output->sval[0];
        mis_config.write_graph = true;
    } else {
        mis_config.write_graph = false;
    }

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}

#endif
