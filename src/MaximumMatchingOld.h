#ifndef MAXIMUM_MATCHING_OLD_H
#define MAXIMUM_MATCHING_OLD_H

#include "SparseArraySet.h"
#include <vector>

class MaximumMatchingOld{
public:
	MaximumMatchingOld();
	MaximumMatchingOld(std::vector<SparseArraySet> &_neighbors);
	void InitialMatching();
	void PPF();
	void MarkReachableVertices();

	std::vector<int> reachableVertices;

protected:
	void MatchAndUpdate(int u);
	int DFS_LA_TS(int startVertex, std::vector<int> &path);
	int DFS_LA_TS_Reverse(int startVertex, std::vector<int> &path);

	// Just for testing
	bool IsMaximalMatching();
	bool CheckUnmatchedLeft();
	bool IsValidVertexCover();
	bool CheckVertexCoverAndMatchingSize();

	std::vector<SparseArraySet> *neighbors;
	std::vector<int> MatchingMates;
	std::vector<int> degree;
	std::vector<int> unmatchedLeft;
	std::vector<int> visited;
	std::vector<int> deg1Vertices;
	std::vector<int> tempUnmatchedVertices;
	std::vector<int> lookAhead;
	std::vector<int> nextNeighbor;
	std::vector<int> nextNeighborLookAhead;
	std::vector<std::vector<int>> augmentingPaths;
	int N;
	int numUnmatchedLeft;
};

#endif //MAXIMUM_MATCHING_H
