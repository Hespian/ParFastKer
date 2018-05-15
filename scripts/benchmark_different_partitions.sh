#!/bin/bash
dirname="../../LinearTimeKernels/"
for filename in $dirname/*.graph; do
    command="../build/benchmark $filename --partition_path=$dirname/partitions/$(basename "$filename")/weight_one_LPA/32.partition --console_log --num_reps=3"
    echo $command
    $command &> ../../results/LinearTimeKernels/lpa/$(basename "$filename")
done
