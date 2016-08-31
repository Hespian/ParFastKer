#include "MaximumMatchingOld.h"
#include <omp.h>
#include <unordered_map>
#include <utility>

MaximumMatchingOld::MaximumMatchingOld(std::vector<SparseArraySet> &_neighbors) : neighbors(&_neighbors), N(_neighbors.size()) {
	MatchingMates = std::vector<int>(2 * N);
	degree = std::vector<int>(N);
	unmatchedLeft = std::vector<int>(2 * N);
	visited = std::vector<int>(2 * N);
	deg1Vertices = std::vector<int>(N);
	reachableVertices = std::vector<int>(2 * N);
	tempUnmatchedVertices = std::vector<int>(2 * N);
	lookAhead = std::vector<int>(2 * N);
	nextNeighbor = std::vector<int>(2 * N);
	nextNeighborLookAhead = std::vector<int> (2 * N);
	long numThreads;
#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	augmentingPaths = std::vector<std::vector<int>>(numThreads);
	for(int i = 0; i < numThreads; i++)
	{
		augmentingPaths[i] = std::vector<int>(2 * N);
	}
}

MaximumMatchingOld::MaximumMatchingOld() {}

void MaximumMatchingOld::PPF() {    
    
#pragma omp parallel for default(shared)
	for(int vertex = 0; vertex < 2 * N; vertex++)
	{
		lookAhead[vertex] = 0;
		nextNeighborLookAhead[vertex] = -1;      
	}
	
	
	int iterations = 0;
    
	while(1)
	{
		iterations++;
		int numTempUnmatchedVertices = 0;
		if(iterations % 2 == 1)
		{
#pragma omp parallel for schedule(static)
			for(int vertex = 0; vertex < 2 * N; vertex++)
			{
				visited[vertex] = 0;
				nextNeighbor[vertex] = -1;
			}
		}
		else
		{
#pragma omp parallel for schedule(static)
			for(int  vertex = 0; vertex < 2 * N; vertex++)
			{
				visited[vertex] = 0;
				nextNeighbor[vertex] = (*neighbors)[vertex].Size();
			}
		}
#pragma omp parallel for schedule(dynamic) default(shared)
		for(int unmatchedVertexIndex = 0; unmatchedVertexIndex < numUnmatchedLeft; unmatchedVertexIndex++)
		{
			
			long threadId = omp_get_thread_num();
			assert(threadId < omp_get_num_threads());
			std::vector<int> &augmentingPath = augmentingPaths[threadId];
			int unmatchedVertex = unmatchedLeft[unmatchedVertexIndex];
			assert(unmatchedVertex < N);
			int augPathLen;
			if(iterations % 2 == 1) // odd iterations
				augPathLen = DFS_LA_TS(unmatchedVertex, augmentingPath) ;
			else
				augPathLen = DFS_LA_TS_Reverse(unmatchedVertex, augmentingPath) ;
			if (augPathLen > 0)
			{
				assert(augPathLen % 2 == 0);
				for(long k=0; k< augPathLen; k += 2)
				{
					MatchingMates[augmentingPath[k]] = augmentingPath[k+1];
					MatchingMates[augmentingPath[k+1]] = augmentingPath[k];
				}
				
			}
			else
			{
				tempUnmatchedVertices[__sync_fetch_and_add(&numTempUnmatchedVertices,1)] = unmatchedVertex;
			}
			
		}
        if( (numTempUnmatchedVertices == 0) || (numUnmatchedLeft == numTempUnmatchedVertices))
        {
        	tempUnmatchedVertices.swap(unmatchedLeft);
			numUnmatchedLeft = numTempUnmatchedVertices;
            break;
        }
		tempUnmatchedVertices.swap(unmatchedLeft);
		numUnmatchedLeft = numTempUnmatchedVertices;
	}
	assert(CheckUnmatchedLeft());
	assert(IsMaximalMatching());
}

int MaximumMatchingOld::DFS_LA_TS(int startVertex, std::vector<int> &path)
{
	
	int top = -1;
	assert(path.size() == 2 * N);
	assert(startVertex >= 0);
	assert(startVertex < N);
	path[++top] = startVertex; // push , path is equivalent to stack
	
	while (top >= 0 )// while stack not empty
	{
		int u = path[top];
		int uDegree = (*neighbors)[u].Size();
		// lookahed part
		while(++nextNeighborLookAhead[u] < uDegree)
		{
			int v = (*neighbors)[u][nextNeighborLookAhead[u]] + N;
			if(__sync_fetch_and_add(&lookAhead[v],1) == 0)
			{
				if(MatchingMates[v] == -1)
				{
					__sync_fetch_and_add(&visited[v],1);
					path[++top] = v; // push
					return top+1; // top = augmenting path length
					
				}
				
			}
		}
		
		while(++nextNeighbor[u] < uDegree)
		{
			
			int v = (*neighbors)[u][nextNeighbor[u]] + N;
			if(__sync_fetch_and_add(&visited[v],1) == 0)
			{
				if(MatchingMates[v] != -1) // means other vertex already allocate this in lookahed phase
				{
					
					path[++top] = v; // push v
					path[++top] = MatchingMates[v]; // push next u
					break;
					
				}
			}
			
		}
		if(nextNeighbor[u] == uDegree)
		{
			top-= 2;// pop
		}
	}
	return top+1;
}





// DFS with lookahead that finds a single augmenting path 
// called from pothen-fan
int MaximumMatchingOld::DFS_LA_TS_Reverse(int startVertex, std::vector<int> &path)
{
	int top = -1;
	assert(path.size() == 2 * N);
	assert(startVertex >= 0);
	assert(startVertex < N);
	path[++top] = startVertex; // push , path is equivalent to stack 
	
	while (top >= 0 )// while stack not empty 
	{
		int u = path[top];
		int uDegree = (*neighbors)[u].Size();
		// lookahed part
		while(++nextNeighborLookAhead[u] < uDegree)
		{
			int v = (*neighbors)[u][nextNeighborLookAhead[u]] + N;
			if(__sync_fetch_and_add(&lookAhead[v],1) == 0)
			{
				if(MatchingMates[v] == -1)
				{
					__sync_fetch_and_add(&visited[v],1);
					path[++top] = v; // push
					return top+1; // top = augmenting path length
					
				}
				
			}
		}
		
		while(--nextNeighbor[u] >= 0)
		{
			
			int v = (*neighbors)[u][nextNeighbor[u]] + N;
			if(__sync_fetch_and_add(&visited[v],1) == 0) 
			{
				if(MatchingMates[v] != -1) // means other vertex already allocate this in lookahed phase 
				{
					
					path[++top] = v; // push v
					path[++top] = MatchingMates[v]; // push next u
					break;
					
				}
			}
			
		}
		if(nextNeighbor[u] == -1)
		{
			top-= 2;// pop
		}
		
		
	}
	return top+1;
}

void MaximumMatchingOld::MarkReachableVertices() {
#pragma omp parallel for
	for(int vertex = 0; vertex < N; ++vertex) {
		reachableVertices[vertex] = 0;
	}
#pragma omp parallel for
	for(int startVertex = 0; startVertex < N; ++startVertex) {
		if(MatchingMates[startVertex] == -1) {
			assert(reachableVertices[startVertex] == 0);
			reachableVertices[startVertex] = 1;
			long threadId = omp_get_thread_num();
			std::vector<int> &stack = augmentingPaths[threadId];
			int top = -1;
			assert(stack.size() == 2 * N);
			assert(startVertex >= 0);
			assert(startVertex < N);
			stack[++top] = startVertex;
			
			while (top >= 0 )
			{
				int vertex = stack[top];
				top --;
				int vertexDegree = (*neighbors)[vertex].Size();
				assert(vertex < N);
				for(int realNeighbor : (*neighbors)[vertex])
				{
					int neighbor = realNeighbor + N;
					int matchingMate = MatchingMates[neighbor];
					assert(matchingMate != -1);
					assert(matchingMate < N);
					if(neighbor == matchingMate)
						continue;

					if(__sync_fetch_and_add(&reachableVertices[neighbor],1) == 0) 
					{
						assert(reachableVertices[matchingMate] == 0);
						reachableVertices[matchingMate] += 1;
						assert(reachableVertices[matchingMate] == 1);
						stack[++top] = matchingMate;
						assert(top < 2 * N);
					}
					
				}				
			}
		}
	}
	assert(IsValidVertexCover());
	assert(CheckVertexCoverAndMatchingSize());
}

// Just for testing
bool MaximumMatchingOld::IsValidVertexCover() {

	for(int vertex = 0; vertex < N; ++vertex) {
		if(reachableVertices[vertex] != 0) {
			for(int realneighbor : (*neighbors)[vertex]) {
				int neighbor = realneighbor + N;
				if(reachableVertices[neighbor] == 0)
					return false;
			}
		}
	}
	return true;
}

// Just for testing
bool MaximumMatchingOld::CheckVertexCoverAndMatchingSize() {
	int vertexCoverSize = 0;
	for(int vertex = 0; vertex < N; ++vertex) {
		if(reachableVertices[vertex] == 0) {
			++vertexCoverSize;
		}
	}

	for(int vertex = N; vertex < 2 * N; ++vertex) {
		if(reachableVertices[vertex] > 0) {
			++vertexCoverSize;
		}
	}

	int matchingSize = 0;
	for(int vertex = 0; vertex < N; ++vertex) {
		if(MatchingMates[vertex] != -1) {
			++matchingSize;
		}
	}

	return vertexCoverSize == matchingSize;
}

void MaximumMatchingOld::MatchAndUpdate(int u)
{
	if(__sync_fetch_and_add(&visited[u],1) != 0) return;	
	for(int v : (*neighbors)[u])
	{
		if(__sync_fetch_and_add(&visited[v + N],1) == 0)
		{
			MatchingMates[u] = v + N;
			MatchingMates[v + N] = u;
			for(int nextNeighbor : (*neighbors)[v])
			{
				if( __sync_fetch_and_add(&degree[nextNeighbor],-1) == 2)
				{					
					MatchAndUpdate(nextNeighbor);
				}		
			}
			break;
		} 
	}
}


void MaximumMatchingOld::InitialMatching()
{	
	numUnmatchedLeft = 0;
	
	double startTime = omp_get_wtime();
	
#pragma omp parallel for default(shared) schedule(static)
	for(int i=0; i < 2*N; i++)
	{
		visited[i] = 0;
		MatchingMates[i] = -1;
	}
	
	int degree1Count = 0;  
	
#pragma omp parallel for default(shared) schedule(static)
	for(int u=0; u < N; u++)
	{      
		degree[u] = (*neighbors)[u].Size();
		if(degree[u] == 1)
		{
			deg1Vertices[__sync_fetch_and_add(&degree1Count,1)] = u;
		}
	}
	
#pragma omp parallel for default(shared)
	for(int u=0; u < degree1Count; u++)
	{
		MatchAndUpdate(deg1Vertices[u]);		  
	}
	
#pragma omp parallel for default(shared) schedule(dynamic,100)
	for(int u=0; u < N; u++)
	{
		if(visited[u] == 0 && degree[u]>0)
			MatchAndUpdate(u);	  
	}
	
#pragma omp parallel for default(shared) 
	for(int u=0; u < N; u++)
	{
		
		if(MatchingMates[u] == -1 && (*neighbors)[u].Size() > 0)
		{
			unmatchedLeft[__sync_fetch_and_add(&numUnmatchedLeft, 1)] = u;
		}
	}

	double endTime = omp_get_wtime();
	double totalTime = endTime - startTime;

	assert(CheckUnmatchedLeft());
	assert(IsMaximalMatching());
}

// Just for testing
bool MaximumMatchingOld::CheckUnmatchedLeft() {
	std::vector<char> markedVertices(N, false);
	for(int unmatchedVertexIndex = 0; unmatchedVertexIndex < numUnmatchedLeft; ++unmatchedVertexIndex) {
		int unmatchedVertex = unmatchedLeft[unmatchedVertexIndex];
		assert(unmatchedVertex < N);
		markedVertices[unmatchedVertex] = true;
	}

	for(int vertex = 0; vertex < N; ++vertex) {
		bool isUnmatched = MatchingMates[vertex] == -1;
		if((*neighbors)[vertex].Size() > 0 && markedVertices[vertex] != isUnmatched)
			return false;
	}
	return true;
}

// Just for testing
bool MaximumMatchingOld::IsMaximalMatching() {
	for(int i = 0; i < numUnmatchedLeft; ++i) {
		int vertex = unmatchedLeft[i];
		for(int neighbor : (*neighbors)[vertex]) {
			if(MatchingMates[neighbor + N] == -1) {
				return false;
			}
		}
	}
	for(int i = 0; i < N; i++) {
		if(MatchingMates[i] != -1) {
			if(MatchingMates[MatchingMates[i]] != i) {
				return false;
			}
		}
	}
	return true;
}
