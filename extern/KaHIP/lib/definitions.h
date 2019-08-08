/******************************************************************************
 * definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef DEFINITIONS_H_CHR
#define DEFINITIONS_H_CHR

#include <limits>
#include <queue>
#include <vector>

#include "limits.h"
#include "macros_assertions.h"
#include "stdio.h"

// allows us to disable most of the output during partitioning
#ifdef KAFFPAOUTPUT
        #define PRINT(x) x
#else
        #define PRINT(x) do {} while (false);
#endif

/**********************************************
 * Constants
 * ********************************************/
//Types needed for the graph ds
#ifdef UNSAFE_LONG
typedef unsigned long 	NodeID;
typedef double 		EdgeRatingType;
typedef unsigned long 	EdgeID;
typedef unsigned long 	PathID;
typedef unsigned long 	PartitionID;
typedef unsigned long 	NodeWeight;
typedef long 		EdgeWeight;
typedef EdgeWeight 	Gain;
typedef long 		Color;
typedef unsigned long 	Count;
typedef std::vector<NodeID> boundary_starting_nodes;
typedef long FlowType;
#else
typedef unsigned int 	NodeID;
typedef double 		EdgeRatingType;
typedef unsigned int 	EdgeID;
typedef unsigned int 	PathID;
typedef unsigned int 	PartitionID;
typedef unsigned int 	NodeWeight;
typedef int 		EdgeWeight;
typedef EdgeWeight 	Gain;
typedef int 		Color;
typedef unsigned int 	Count;
typedef std::vector<NodeID> boundary_starting_nodes;
typedef long FlowType;
#endif // UNSAFE_LONG

const EdgeID UNDEFINED_EDGE            = std::numeric_limits<EdgeID>::max();
const NodeID NOTMAPPED                 = std::numeric_limits<EdgeID>::max();
const NodeID UNDEFINED_NODE            = std::numeric_limits<NodeID>::max();
const PartitionID INVALID_PARTITION    = std::numeric_limits<PartitionID>::max();
const PartitionID BOUNDARY_STRIPE_NODE = std::numeric_limits<PartitionID>::max();
const int NOTINQUEUE 		       = std::numeric_limits<int>::max();
const int ROOT 			       = 0;

//for the gpa algorithm
struct edge_source_pair {
        EdgeID e;
        NodeID source;       
};

struct source_target_pair {
        NodeID source;       
        NodeID target;       
};

//matching array has size (no_of_nodes), so for entry in this table we get the matched neighbor
typedef std::vector<NodeID> CoarseMapping;
typedef std::vector<NodeID> Matching;
typedef std::vector<NodeID> NodePermutationMap;

typedef double ImbalanceType;
//Coarsening
typedef enum {
        EXPANSIONSTAR, 
        EXPANSIONSTAR2, 
 	WEIGHT, 
	PSEUDOGEOM, 
	EXPANSIONSTAR2ALGDIST, 
} EdgeRating;

typedef enum {
        PERMUTATION_QUALITY_NONE, 
	PERMUTATION_QUALITY_FAST,  
	PERMUTATION_QUALITY_GOOD
} PermutationQuality;

typedef enum {
        MATCHING_RANDOM, 
	MATCHING_GPA, 
	MATCHING_RANDOM_GPA,
        CLUSTER_COARSENING
} MatchingType;

typedef enum {
	INITIAL_PARTITIONING_RECPARTITION, 
	INITIAL_PARTITIONING_BIPARTITION
} InitialPartitioningType;

typedef enum {
        REFINEMENT_SCHEDULING_FAST, 
	REFINEMENT_SCHEDULING_ACTIVE_BLOCKS, 
	REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY
} RefinementSchedulingAlgorithm;

typedef enum {
        REFINEMENT_TYPE_FM, 
	REFINEMENT_TYPE_FM_FLOW, 
	REFINEMENT_TYPE_FLOW
} RefinementType;

typedef enum {
        STOP_RULE_SIMPLE, 
	STOP_RULE_MULTIPLE_K, 
	STOP_RULE_STRONG 
} StopRule;

typedef enum {
        BIPARTITION_BFS, 
	BIPARTITION_FM
} BipartitionAlgorithm ;

typedef enum {
        KWAY_SIMPLE_STOP_RULE, 
	KWAY_ADAPTIVE_STOP_RULE
} KWayStopRule;

typedef enum {
        COIN_RNDTIE, 
	COIN_DIFFTIE, 
	NOCOIN_RNDTIE, 
	NOCOIN_DIFFTIE 
} MLSRule;

typedef enum {
        CYCLE_REFINEMENT_ALGORITHM_PLAYFIELD, 
        CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL, 
	CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL_PLUS
} CycleRefinementAlgorithm;

typedef enum {
        RANDOM_NODEORDERING, 
        DEGREE_NODEORDERING
} NodeOrderingType;


#endif

