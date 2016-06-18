#!/bin/bash

cd ../build
graph_base_dir="/global_data/c_schulz/graph_collection/all_dimacs/dimacs10/archive/data/streets/"
graph="europe.osm.graph"
sizes=("1" "2" "4" "8" "16" "32" "48" "64")
modes=("kahip" "parallel_kahip" "lpa")

cmake .. -DNO_PREPROCESSING=OFF
make -j 4

for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode-preprocessing
	done
done

cmake .. -DNO_PREPROCESSING=ON
make -j 4
for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode
	done
done

graph_base_dir="/global_data/c_schulz/graph_collection/rhgs_large/"
graph="RHG-100000000-nodes-1000000000-edges.graph"

cmake .. -DNO_PREPROCESSING=OFF
make -j 4

for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode-preprocessing
	done
done

cmake .. -DNO_PREPROCESSING=ON
make -j 4
for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode
	done
done

graph="RHG-100000000-nodes-2000000000-edges.graph"

cmake .. -DNO_PREPROCESSING=OFF
make -j 4

for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode-preprocessing
	done
done

cmake .. -DNO_PREPROCESSING=ON
make -j 4
for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode
	done
done

graph_base_dir="/global_data/c_schulz/graph_collection/mis_huge/"
graph="it-2004-sorted.graph"

cmake .. -DNO_PREPROCESSING=OFF
make -j 4

for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode-preprocessing
	done
done

cmake .. -DNO_PREPROCESSING=ON
make -j 4
for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode
	done
done

graph="sk-2005-sorted.graph"

cmake .. -DNO_PREPROCESSING=OFF
make -j 4

for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode-preprocessing
	done
done

cmake .. -DNO_PREPROCESSING=ON
make -j 4
for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode &> ../../results/$graph-$size-$mode
	done
done
