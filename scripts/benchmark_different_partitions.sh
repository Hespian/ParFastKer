#!/bin/bash
dirname="../../graphs"
for filename in $dirname/*.graph; do
    command="../build/benchmark $filename --partition_path=$dirname/partitions/$(basename "$filename")/ --console_log --num_reps=1"
    echo $command
    $command &> ../../results/convergence/LP_LAST_UNCONFINED_6/$(basename "$filename")
done
