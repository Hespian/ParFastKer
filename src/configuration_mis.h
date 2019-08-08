/**
 * configuration.h
 * Purpose: Contains preset configurations for the algorithm.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
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

#ifndef _CONFIGURATION_MIS_H_
#define _CONFIGURATION_MIS_H_

#include "definitions.h"
#include "mis_config.h"

class configuration_mis {
    public:
        /**
         * Default Constructor.
         */
        configuration_mis() {} ;

        /**
         * Default Destructor.
         */
        virtual ~configuration_mis() {};

        /**
         * Set the standard configuration.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void standard( MISConfig & config );

        /**
         * Set the configuration for social network graphs.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void social( MISConfig & config );

        /**
         * Set the configuration for the experimental evaluation.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void full_standard( MISConfig & config );

        /**
         * Set the configuration for the experimental evaluation 
         * for social network graphs.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void full_social( MISConfig & config );
};

inline void configuration_mis::standard( MISConfig & mis_config ) {
    // Reductions
    mis_config.all_reductions                         = true;
    mis_config.num_reps                               = 5;
}

inline void configuration_mis::social( MISConfig & mis_config ) {
    standard(mis_config);
}

inline void configuration_mis::full_standard( MISConfig & mis_config ) {
    standard(mis_config);
}

inline void configuration_mis::full_social( MISConfig & mis_config ) {
    full_standard(mis_config);
}

#endif 
