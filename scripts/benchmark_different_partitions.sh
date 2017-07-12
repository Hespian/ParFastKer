#!/bin/bash
dirname="../../graphs"
for filename in $dirname/*.graph; do
    command="../build/benchmark $filename --partition_path=$dirname/partitions/$(basename "$filename")/ --console_log --num_reps=1"
    echo $command
    $command &> ../../results/convergence/SMALL_SLOPE_GLOBAL/$(basename "$filename")
done
