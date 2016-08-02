import matplotlib.pyplot as plt
from collections import namedtuple
import numpy as np
import sys
import os
import seaborn as sns

inputFile = sys.argv[1]
plotsDir = inputFile + "-distribution.png"

weights = []
file = open(inputFile, "r")
for line in file:
	weight = int(line)
	weights.append(weight)


plt.hist(weights, bins = 20)
plt.title('distribution of weights')
plt.yscale('log')
plt.savefig(plotsDir)