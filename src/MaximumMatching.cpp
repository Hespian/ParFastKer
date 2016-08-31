#include <parallel/numeric>
#include "MaximumMatching.h"

MaximumMatching::MaximumMatching(std::vector<std::vector<int>> const &adjacencyArray) {
	g = (graph *) malloc(sizeof(graph));
	g->weight = NULL;
	long numVertices = adjacencyArray.size();
	g->vtx_pointer = new long[numVertices * 2 + 1];
	degrees = std::vector<long>(numVertices * 2);
	#pragma omp parallel for
	for(int i = 0; i < numVertices; ++i ) {
		degrees[i] = adjacencyArray[i].size();
		degrees[i + numVertices] = adjacencyArray[i].size();
	}
	auto end_ptr = __gnu_parallel::partial_sum(degrees.begin(), degrees.end(), g->vtx_pointer);
	assert(end_ptr == &(g->vtx_pointer[2 * numVertices]));
	long numEdges = g->vtx_pointer[2 * numVertices - 1];
	#pragma omp parallel for
	for(int i = 0; i < 2 * numVertices; ++i ) {
		g->vtx_pointer[i] -= degrees[i];
	}
	g->vtx_pointer[2 * numVertices] = numEdges;
	assert(g->vtx_pointer[2 * numVertices] == g->vtx_pointer[2 * numVertices - 1] + degrees[2 * numVertices - 1]);
	g->n = 2 * numVertices;
	g->m = numEdges;
	g->nrows = numVertices;
	g->endV = new long[numEdges];
	#pragma omp parallel for
	for(int i = 0; i < numVertices; ++i) {
		long index = g->vtx_pointer[i];
		for(int neighbor: adjacencyArray[i]) {
			g->endV[index++] = neighbor + numVertices;
		}
		index = g->vtx_pointer[i + numVertices];
		for(int neighbor: adjacencyArray[i]) {
			g->endV[index++] = neighbor;
		}
	}
}

void MaximumMatching::LoadGraph(std::vector<SparseArraySet> &neighbors) {
	assert(neighbors.size() == g->nrows);
	#pragma omp parallel for
	for(int i = 0; i < g->nrows; ++i ) {
		degrees[i] = neighbors[i].Size();
		degrees[i + g->nrows] = neighbors[i].Size();
	}
	auto end_ptr = __gnu_parallel::partial_sum(degrees.begin(), degrees.end(), g->vtx_pointer);
	assert(end_ptr == &(g->vtx_pointer[g->n]));
	long numEdges = g->vtx_pointer[g->n - 1];
	#pragma omp parallel for
	for(int i = 0; i < g->n; ++i ) {
		g->vtx_pointer[i] -= degrees[i];
	}
	g->vtx_pointer[g->n] = numEdges;
	assert(g->vtx_pointer[g->n] == g->vtx_pointer[g->n - 1] + degrees[g->n - 1]);
	g->m = numEdges;
	
	#pragma omp parallel for
	for(int i = 0; i < g->nrows; ++i) {
		long index = g->vtx_pointer[i];
		for(int neighbor: neighbors[i]) {
			g->endV[index++] = neighbor + g->nrows;
		}
		index = g->vtx_pointer[i + g->nrows];
		for(int neighbor: neighbors[i]) {
			g->endV[index++] = neighbor;
		}
	}
}

MaximumMatching::~MaximumMatching() {
	free_graph(g);
	free(g);
}

void free_graph( graph* bGraph)
{
	delete [] bGraph->vtx_pointer;
	delete [] bGraph->endV;
    if(bGraph->weight) delete [] bGraph->weight;
}
