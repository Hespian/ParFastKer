#!/bin/bash
for filename in ../../graphs/*.graph; do
    command="../build/benchmark $filename --partition_path=../../graphs/partitions/$(basename "$filename")/ --console_log --num_reps=1"
    echo $command
    $command &> ../../no_diamond_no_twin_unlimited_isolated_cliques/$(basename "$filename")
done
