#ifndef PROFILING_HELPER_H
#define PROFILING_HELPER_H

#include "SparseArraySet.h"

#include <time.h>
#include <unordered_map>

#ifdef PROFILING
struct ProfilingHelper_t {
	std::vector<clock_t> startTimersPerPartition;
	std::vector<int> currentVertexPartition;
	std::vector<clock_t> startTimersUpdateNeighborhoodPerPartition;
	std::vector<int> currentVertexUpdateNeighborhoodPartition;

	std::vector<SparseArraySet> *neighborsPtr;

	std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUpdateNeighborhoodsPerDegree; // partition, degree

    std::vector<std::vector<clock_t>> timesUnsuccessfulFoldDegree; //partition
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulFoldAdjacentPerNeighborDegree; // partition, first neighbor degree
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize; // partition, size of 2-neighborhood
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesSuccessfulFoldPerTwoNeighborhoodSize; // partition, size of 2-neighborhood

    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree; // partition, degree
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulIsolatedCliqueNoCliquePerDegree; // partition, degree
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesSuccessfulIsolatedCliquePerDegree; // partition, degree
};

namespace {

	int twoNeighbordhoodSize(ProfilingHelper_t *profilingHelper, int const vertex) {
		int size = (*profilingHelper->neighborsPtr)[vertex].Size();
		for(int neighbor : (*profilingHelper->neighborsPtr)[vertex]) {
	        size += (*profilingHelper->neighborsPtr)[neighbor].Size();
	    }
	    return size;
	}

	void addNewTime(std::unordered_map<int, std::vector<clock_t>> &map, int const index, clock_t const time) {
		map[index].push_back(time);
	}

	void printVectorMapVector(std::vector<std::unordered_map<int, std::vector<clock_t>>> &vec, const char name[]) {
		std::cout << "#########################################################" << std::endl;
		for(int partition = 0; partition < vec.size(); partition++) {
			std::cout << "---------------------------------------------------------" << std::endl;
			for(std::pair<const int, std::vector<clock_t>> element: vec[partition]) {
				int size = element.first;
				std::vector<clock_t> times = element.second;
				std::cout << name << "[" <<partition << "]: " << size << " -";
				for(clock_t time : times) {
					std::cout << " " << time;
				}
				std::cout << std::endl;
			}
		}
	}

	void printclock_tVector(std::vector<std::vector<clock_t>> &vec, const char name[]) {
		std::cout << "#########################################################" << std::endl;
		for(int partition = 0; partition < vec.size(); partition++) {
			std::cout << name << "[" <<partition << "]:";
			for(clock_t time : vec[partition]) {
				std::cout << " " << time;
			}
			std::cout << std::endl;
		}
	}
}

void profilingInit(ProfilingHelper_t *profilingHelper, std::vector<SparseArraySet> *neighbors, int const numPartitions)  {
	profilingHelper->neighborsPtr = neighbors;
	profilingHelper->timesUpdateNeighborhoodsPerDegree = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->timesUnsuccessfulFoldDegree = std::vector<std::vector<clock_t>>(numPartitions);
	profilingHelper->timesUnsuccessfulFoldAdjacentPerNeighborDegree = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->timesSuccessfulFoldPerTwoNeighborhoodSize = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->timesUnsuccessfulIsolatedCliqueNoCliquePerDegree = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->timesSuccessfulIsolatedCliquePerDegree = std::vector<std::unordered_map<int, std::vector<clock_t>>>(numPartitions);
	profilingHelper->startTimersPerPartition = std::vector<clock_t>(numPartitions);
	profilingHelper->currentVertexPartition = std::vector<int>(numPartitions);
	profilingHelper->startTimersUpdateNeighborhoodPerPartition = std::vector<clock_t>(numPartitions);
	profilingHelper->currentVertexUpdateNeighborhoodPartition = std::vector<int>(numPartitions);
}

void profilingStartClock(ProfilingHelper_t *profilingHelper, int const partition, int const vertex) {
	profilingHelper->currentVertexPartition[partition] = vertex;
	profilingHelper->startTimersPerPartition[partition] = clock();
}

void profilingStartClockUpdateNeighborhood(ProfilingHelper_t *profilingHelper, int const partition, int const vertex) {
	profilingHelper->currentVertexUpdateNeighborhoodPartition[partition] = vertex;
	profilingHelper->startTimersUpdateNeighborhoodPerPartition[partition] = clock();
}

void profilingAddTimeUpdateNeighborhood(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersUpdateNeighborhoodPerPartition[partition];
	int vertex = profilingHelper->currentVertexUpdateNeighborhoodPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesUpdateNeighborhoodsPerDegree[partition], degree, time);
}

void profilingAddTimeUnsuccessfulFoldDegree(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	profilingHelper->timesUnsuccessfulFoldDegree[partition].push_back(time);
}

void profilingAddTimeUnsuccessfulFoldAdjacent(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	int vertex = profilingHelper->currentVertexPartition[partition];
	int neighbor = (*profilingHelper->neighborsPtr)[vertex][0];
	int neighborDegree = (*profilingHelper->neighborsPtr)[neighbor].Size();
	addNewTime(profilingHelper->timesUnsuccessfulFoldAdjacentPerNeighborDegree[partition], neighborDegree, time);
}

void profilingAddTimeUnsuccessfulFoldWrongPartition(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	int vertex = profilingHelper->currentVertexPartition[partition];
	int NeighbordhoodSize = twoNeighbordhoodSize(profilingHelper, vertex);
	addNewTime(profilingHelper->timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize[partition], NeighbordhoodSize, time);
}

void profilingAddTimeSuccessfulFold(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	int vertex = profilingHelper->currentVertexPartition[partition];
	int NeighbordhoodSize = twoNeighbordhoodSize(profilingHelper, vertex);
	addNewTime(profilingHelper->timesSuccessfulFoldPerTwoNeighborhoodSize[partition], NeighbordhoodSize, time);
}

void profilingAddTimeUnsuccessfulIsolatedCliqueDegreeOrPartition(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	int vertex = profilingHelper->currentVertexPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree[partition], degree, time);
}

void profilingAddTimeUnsuccessfulIsolatedCliqueNoClique(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	int vertex = profilingHelper->currentVertexPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesUnsuccessfulIsolatedCliqueNoCliquePerDegree[partition], degree, time);
}

void profilingAddTimeSuccessfulIsolatedClique(ProfilingHelper_t *profilingHelper, int const partition) {
	clock_t time = clock() - profilingHelper->startTimersPerPartition[partition];
	int vertex = profilingHelper->currentVertexPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesSuccessfulIsolatedCliquePerDegree[partition], degree, time);
}

void profilingPrint(ProfilingHelper_t *profilingHelper) {
	std::cout << "#########################################################" << std::endl;
	std::cout << std::endl << "Profiling:" << std::endl;

	printVectorMapVector(profilingHelper->timesUpdateNeighborhoodsPerDegree, "UpdateNeighborhoodsPerDegree");

	printclock_tVector(profilingHelper->timesUnsuccessfulFoldDegree, "UnsuccessfulFoldDegree");
	printVectorMapVector(profilingHelper->timesUnsuccessfulFoldAdjacentPerNeighborDegree, "UnsuccessfulFoldAdjacentPerNeighborDegree");
	printVectorMapVector(profilingHelper->timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize, "UnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize");
	printVectorMapVector(profilingHelper->timesSuccessfulFoldPerTwoNeighborhoodSize, "SuccessfulFoldPerTwoNeighborhoodSize");

	printVectorMapVector(profilingHelper->timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree, "UnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree");
	printVectorMapVector(profilingHelper->timesUnsuccessfulIsolatedCliqueNoCliquePerDegree, "UnsuccessfulIsolatedCliqueNoCliquePerDegree");
	printVectorMapVector(profilingHelper->timesSuccessfulIsolatedCliquePerDegree, "SuccessfulIsolatedCliquePerDegree");

	std::cout << "#########################################################" << std::endl;
}
#else
struct ProfilingHelper_t {
	char nothing[0];
};
#define profilingInit (void)sizeof
#define profilingStartClock (void)sizeof
#define profilingStartClockUpdateNeighborhood (void)sizeof
#define profilingAddTimeUpdateNeighborhood (void)sizeof
#define profilingAddTimeUnsuccessfulFoldDegree (void)sizeof
#define profilingAddTimeUnsuccessfulFoldAdjacent (void)sizeof
#define profilingAddTimeUnsuccessfulFoldWrongPartition (void)sizeof
#define profilingAddTimeSuccessfulFold (void)sizeof
#define profilingAddTimeUnsuccessfulIsolatedCliqueDegreeOrPartition (void)sizeof
#define profilingAddTimeUnsuccessfulIsolatedCliqueNoClique (void)sizeof
#define profilingAddTimeSuccessfulIsolatedClique (void)sizeof
#define profilingPrint (void)sizeof
#endif
#endif //PROFILING_HELPER_H
