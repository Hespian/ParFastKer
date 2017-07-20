#!/bin/bash
dirname="../../../triangle_counting_paper/MIS_sigmod_pub/results/LinearTimeKernels"
for filename in $dirname/*.graph; do
    command="../build/benchmark $filename --partition_path=$dirname/partitions/$(basename "$filename")/ --console_log --num_reps=3"
    echo $command
    $command &> ../../results/LinearTimeKernels/$(basename "$filename")
done
