# ParFastKer

This is the Code written for the paper "Scalable Kernelization for Maximum Independent Sets" if you would like to use it for your publication, please cite the following paper:

Hespe, Demian, Christian Schulz, and Darren Strash. "Scalable kernelization for maximum independent sets." Journal of Experimental Algorithmics (JEA) 24.1 (2019): 1-16.

In the above paper, we first kernelize the graph using the LinearTime algorithm from https://github.com/LijunChang/Near-Maximum-Independent-Set
 
For partitioning, we used ParHIP (https://github.com/schulzchristian/KaHIP)

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

## License

All files are under the MIT license, except for  src/MaximumMatching.cpp which was released under the BSD 3-clause license by the original authors
