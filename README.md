# ParFastKer

This is the Code written for the paper "Scalable Kernelization for Maximum Independent Sets" if you would like to use it for your publication, please cite the following paper:

Hespe, Demian, Christian Schulz, and Darren Strash. "Scalable kernelization for maximum independent sets." 2018 Proceedings of the Twentieth Workshop on Algorithm Engineering and Experiments (ALENEX). Society for Industrial and Applied Mathematics, 2018.

It was not written for production use but if you want to use it to reduce a graph quickly or run benchmarks for your own, the following instructions should help with that.

## Installation


```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## Usage

```bash
./benchmark [input file] --partition_path=[partition file] --output=[output file] --console_log
```

[input file] should be an unweighted graph in the METIS graph format

[partition file] should be a file containing one line for each vertex in the graph specifying it's partition index, with the first index being 0. The file name should have the form [number of blocks].partition
