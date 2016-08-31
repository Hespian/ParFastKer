#ifndef MAXIMUM_MATCHING_H
#define MAXIMUM_MATCHING_H

#include "SparseArraySet.h"
#include <vector>

typedef struct /* the bipartite graph data structure */
{
	long n; // numver of vertices in both sides
	long nrows; // number of vertices in the left side
	long m; // number of edges
	long* vtx_pointer; // an array of size n+1 storing the pointer in endV array
    long* endV; //an array of size m that stores the second vertex of an edge.
	double* weight; // not used in unweighted graph
} graph;

void free_graph (graph* bGraph);

class MaximumMatching{
public:

	MaximumMatching(std::vector<std::vector<int>> const &adjacencyArray);
	~MaximumMatching();
	void LoadGraph(std::vector<SparseArraySet> &neighbors);

	std::vector<int> reachableVertices;

protected:
	graph* g;
	std::vector<long> degrees;

};

#endif //MAXIMUM_MATCHING_H
