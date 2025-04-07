#-------------------------------------------------------------------------------
#                               IMPORTS
#-------------------------------------------------------------------------------

import os
import ast
import time
import numpy as np
from matplotlib import pyplot as plt

#-------------------------------------------------------------------------------
#                               CONSTANTS
#-------------------------------------------------------------------------------

MIN_N = 2
MAX_N = 50
SHOTS = 512
REPETITIONS = 10
CIRCUIT_SHOT_TIME = 4.2
ANNEALING_SHOT_TIME = 0.0001
DIRECTORY = '../data/circuits/'
DOMAIN = np.arange(MIN_N, MAX_N+1)
CLASSICAL_TIME = [0] * len(DOMAIN)
CIRCUIT_TIME = [CIRCUIT_SHOT_TIME] * (len(DOMAIN))
ANNEALING_TIME = [ANNEALING_SHOT_TIME * SHOTS] * (len(DOMAIN))

with open('solutionTimes.txt', 'w+') as file:
    file.write(f'Annealing time: {ANNEALING_TIME}\n')

#-------------------------------------------------------------------------------
#                               FUNCTIONS
#-------------------------------------------------------------------------------

def f(x):
    if x[0] == '0':
        return x
    else:
        return ''.join(map(str, np.array([1] * len(x)) - np.array(list(map(int, list(x))))))

#-------------------------------------------------------------------------------
#                           CLASSICAL ALGORITHM
#-------------------------------------------------------------------------------

print('Classical')
for n in DOMAIN:
    print(n)
    index = n - MIN_N
    startTime = time.time()
    successes = 0
    while successes < REPETITIONS:
        evals = dict()
        done = False
        while not done:
            x = ''.join(np.random.choice(['0', '1'], size=n, replace=True))
            fx = f(x)
            if fx in evals.keys():
                if evals[fx] != x:
                    done = True
            else:
                evals[fx] = x
        successes += 1
    CLASSICAL_TIME[index] = (time.time() - startTime) / REPETITIONS
    print(time.time() - startTime)

with open('solutionTimes.txt', 'a') as file:
    file.write(f'Classical time: {CLASSICAL_TIME}\n')

#-------------------------------------------------------------------------------
#                           CIRCUIT ALGORITHM
#-------------------------------------------------------------------------------

print('Circuit')
for filename in os.listdir(DIRECTORY):
    with open(DIRECTORY + filename) as file:
        startTime = time.time()
        successes = 0
        n = int(filename.split('-')[1][:-4])
        print(n)
        index = n - MIN_N
        counts = ast.literal_eval(file.readlines()[0])
        try:
            del counts['0'*n]
        except:
            pass
        probs = np.array(list(counts.values())) / sum(list(counts.values()))
        while successes < REPETITIONS:
            subsample = np.random.choice(list(counts.keys()), size=n-1, replace=False, p=probs)
            subsystem = np.array([[int(i) for i in result] for result in subsample])
            if np.all(np.sum(subsystem, axis=0) >= np.ones(n)) and np.all(np.sum(subsystem, axis=1)%2 == np.zeros(n-1)):
                successes += 1
                print((time.time() - startTime))
        CIRCUIT_TIME[index] += ((time.time() - startTime) / REPETITIONS)

with open('solutionTimes.txt', 'a') as file:
    file.write(f'Circuit time: {CIRCUIT_TIME}\n')

#-------------------------------------------------------------------------------
#                               PLOTTING
#-------------------------------------------------------------------------------

plt.figure(figsize=(10,10))
plt.semilogy(DOMAIN, ANNEALING_TIME, label='Annealing')
plt.semilogy(DOMAIN, CIRCUIT_TIME, label='Circuit')
plt.semilogy(DOMAIN, CLASSICAL_TIME, label='Classical')
plt.legend(loc='best')
plt.xlabel('$n$')
plt.ylabel('Time (s)')
plt.savefig('timeComparison.png', dpi=300)