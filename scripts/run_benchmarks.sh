#!/bin/bash

cd fast-reductions/build
graph_base_dir="../../../graphs/from_darren/"
graph="petster-cat.graph"
sizes=("1" "2")
modes=("kahip" "parallel_kahip" "lpa")

cmake .. -DNO_PREPROCESSING=OFF
make -j 4

for size in "${sizes[@]}"
do
   :
	for mode in "${modes[@]}"
	do
	   :
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode > ../../results/$graph-$size-$mode-preprocessing
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
		./benchmark $graph_base_dir$graph --console_log --config social --num_partitions=$size --partitioner=$mode > ../../results/$graph-$size-$mode
	done
done