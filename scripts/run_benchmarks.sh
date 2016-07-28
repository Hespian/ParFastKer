#!/bin/bash
python partition_graphs.py ../../graphs > partition_script.sh
python benchmark_directory.py ../../graphs/ 5 > benchmark_script.sh
chmod +x partition_script.sh
./partition_script.sh
chmod +x benchmark_script.sh
./benchmark_script.sh