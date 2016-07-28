import sys
import os
import subprocess

graphDir = sys.argv[1]
numPartitionSet = [1, 2, 4, 8, 16, 32, 48, 64]

def makedir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

for file in os.listdir(graphDir):
    if file.endswith(".graph"):
        graphPath = os.path.join(graphDir, file)
        partitionsDir = os.path.join(graphDir, "partitions", file)
        makedir(partitionsDir)
        targetDir = os.path.join(partitionsDir, "weight_one")
        makedir(targetDir)
        for numPartitions in numPartitionSet:
            targetFile = os.path.join(targetDir, str(numPartitions) + ".partition")
            if not os.path.exists(targetFile):
                outputDumpFile = os.path.join(targetDir, "partition_output")
                callString = "mpirun -n 16 ../../parallel_social_partitioning_package_weighted/deploy/parallel_label_compress " + graphPath + " --k=" + str(numPartitions) + " --preconfiguration=ultrafast --seed 1337 > " + outputDumpFile
                print(callString)
                subprocess.call(callString, shell=True)
                outputFile = os.path.join(os.getcwd(), "tmppartition")
                print("Renaming " + outputFile + " to " + targetFile)
                os.rename(outputFile, targetFile)
