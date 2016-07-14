#ifndef PROFILING_HELPER_H
#define PROFILING_HELPER_H

#include "SparseArraySet.h"

#include <time.h>
#include <unordered_map>

class ProfilingHelper {
public:
	ProfilingHelper(std::vector<SparseArraySet> &neighbors, int const numPartitions) 
	: neighborsRef(&neighbors) 
	, timesUpdateNeighborhoodsPerDegree(numPartitions)
	, timesUnsuccessfulFoldDegree(numPartitions)
	, timesUnsuccessfulFoldAdjacentPerNeighborDegree(numPartitions)
	, timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize(numPartitions)
	, timesSuccessfulFoldPerTwoNeighborhoodSize(numPartitions)
	, timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree(numPartitions)
	, timesUnsuccessfulIsolatedCliqueNoCliquePerDegree(numPartitions)
	, timesSuccessfulIsolatedCliquePerDegree(numPartitions)
	{}

	ProfilingHelper() {}

	void addTimeUpdateNeighborhood(int const partition, int const vertex, clock_t const time) {
		int degree = (*neighborsRef)[vertex].Size();
		addNewTime(timesUpdateNeighborhoodsPerDegree[partition], degree, time);
	}

	void addTimeUnsuccessfulFoldDegree(int const partition, clock_t const time) {
		timesUnsuccessfulFoldDegree[partition].push_back(time);
	}

	void addTimeUnsuccessfulFoldAdjacent(int const partition, int const vertex, clock_t const time) {
		int neighbor = (*neighborsRef)[vertex][0];
		int neighborDegree = (*neighborsRef)[neighbor].Size();
		addNewTime(timesUnsuccessfulFoldAdjacentPerNeighborDegree[partition], neighborDegree, time);
	}

	void addTimeUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize(int const partition, int const vertex, clock_t const time) {
		int NeighbordhoodSize = twoNeighbordhoodSize(vertex);
		addNewTime(timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize[partition], NeighbordhoodSize, time);
	}

	void addTimeSuccessfulFoldPerTwoNeighborhoodSize(int const partition, int const vertex, clock_t const time) {
		int NeighbordhoodSize = twoNeighbordhoodSize(vertex);
		addNewTime(timesSuccessfulFoldPerTwoNeighborhoodSize[partition], NeighbordhoodSize, time);
	}

	void addTimeUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree(int const partition, int const vertex, clock_t const time) {
		int degree = (*neighborsRef)[vertex].Size();
		addNewTime(timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree[partition], degree, time);
	}

	void addTimeUnsuccessfulIsolatedCliqueNoCliquePerDegree(int const partition, int const vertex, clock_t const time) {
		int degree = (*neighborsRef)[vertex].Size();
		addNewTime(timesUnsuccessfulIsolatedCliqueNoCliquePerDegree[partition], degree, time);
	}

	void addTimeSuccessfulIsolatedCliquePerDegree(int const partition, int const vertex, clock_t const time) {
		int degree = (*neighborsRef)[vertex].Size();
		addNewTime(timesSuccessfulIsolatedCliquePerDegree[partition], degree, time);
	}

	void print() {
		std::cout << "#########################################################" << std::endl;
		std::cout << std::endl << "Profiling:" << std::endl;

		printVectorMapVector(timesUpdateNeighborhoodsPerDegree, "UpdateNeighborhoodsPerDegree");

		printclock_tVector(timesUnsuccessfulFoldDegree, "UnsuccessfulFoldDegree");
		printVectorMapVector(timesUnsuccessfulFoldAdjacentPerNeighborDegree, "UnsuccessfulFoldAdjacentPerNeighborDegree");
		printVectorMapVector(timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize, "UnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize");
		printVectorMapVector(timesSuccessfulFoldPerTwoNeighborhoodSize, "SuccessfulFoldPerTwoNeighborhoodSize");

		printVectorMapVector(timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree, "UnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree");
		printVectorMapVector(timesUnsuccessfulIsolatedCliqueNoCliquePerDegree, "UnsuccessfulIsolatedCliqueNoCliquePerDegree");
		printVectorMapVector(timesSuccessfulIsolatedCliquePerDegree, "SuccessfulIsolatedCliquePerDegree");

		std::cout << "#########################################################" << std::endl;
	}

protected:

	int twoNeighbordhoodSize(int const vertex) {
		int size = (*neighborsRef)[vertex].Size();
		for(int neighbor : (*neighborsRef)[vertex]) {
	        size += (*neighborsRef)[neighbor].Size();
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

	std::vector<SparseArraySet> *neighborsRef;

	std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUpdateNeighborhoodsPerDegree; // partition, degree

    std::vector<std::vector<clock_t>> timesUnsuccessfulFoldDegree; //partition
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulFoldAdjacentPerNeighborDegree; // partition, first neighbor degree
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulFoldWrongPartitionPerTwoNeighborhoodSize; // partition, size of 2-neighborhood
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesSuccessfulFoldPerTwoNeighborhoodSize; // partition, size of 2-neighborhood

    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulIsolatedCliqueDegreeOrPartitionPerDegree; // partition, degree
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesUnsuccessfulIsolatedCliqueNoCliquePerDegree; // partition, degree
    std::vector<std::unordered_map<int, std::vector<clock_t>>> timesSuccessfulIsolatedCliquePerDegree; // partition, degree
};
#endif //PROFILING_HELPER_H
