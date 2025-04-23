import os
import ast
import numpy as np
from matplotlib import pyplot as plt

DIRECTORY = '../data/circuits/'
DOMAIN = np.arange(2, 51, 1)
GATE_COUNTS = np.zeros(49)

for filename in os.listdir(DIRECTORY):
    with open(DIRECTORY + filename) as file:
        index = int(filename.split('-')[1][:-4]) - 2
        gateDict = ast.literal_eval(file.readlines()[2][12:-1])
        GATE_COUNTS[index] = sum(gateDict.values())

plt.figure(figsize=(10,10))
plt.plot(DOMAIN, GATE_COUNTS)
plt.xlabel('$n$')
plt.ylabel('Gate Count')
plt.savefig('gateCounts.png', dpi=300)