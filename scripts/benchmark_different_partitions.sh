#!/bin/bash
for filename in ../../graphs/*.graph; do
    command="../build/benchmark $filename --partition_path=../../graphs/partitions/$(basename "$filename")/ --console_log --num_reps=3"
    echo $command
    $command &> ../../results/overpartitioning_no_LP_max_32/$(basename "$filename")
done
