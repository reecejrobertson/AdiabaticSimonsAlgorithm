import dimod
import matplotlib.pyplot as plt
from collections import Counter
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

# TODO: Generalize this to larger problem [goal: n=20]

#TODO: For the moderately sized working problem:
    # TODO: 3d plot with annealing times [10ms, 20ms, 50ms, 100ms, 200ms, 500ms, 1000ms]
    # TODO: Try this from all starting states for reverse annealing
        # TODO: Add pausing into schedule
    # TODO: Repeat with forward annealing & fast annealing
    # TODO: 8192 samples up to 200ms, 1024 samples for 500-1000ms
    # TODO: Vary s3 variable

# # make QUBO for the Pm problem

# from copy import deepcopy
# import itertools
# from operator import itemgetter
# import numpy as np

# import matplotlib
# import matplotlib.pyplot as plt

# import neal
# import dimod

# # these are D-Wave modules

# from dwave.system import EmbeddingComposite, DWaveSampler, LeapHybridSampler
# from dwave.system.composites import FixedEmbeddingComposite
# from minorminer import find_embedding

# # Problem input

# s = neal.SimulatedAnnealingSampler()
# sampleset = s.sample_qubo(
#     Q, beta_range = (0.01, 10), num_sweeps = 200,
#     num_reads = no_runs, beta_schedule_type="geometric"
# )

# https://github.com/iitis/parallel_machines/


# QPU in Juelich
# qpu_sampler = DWaveSampler(solver='Advantage_system5.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4',region="eu-central-1")  # Pegasus Germany
# qpu_sampler = DWaveSampler(solver='Advantage_system4.1',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Pegasus 
# qpu_sampler = DWaveSampler(solver='Advantage_system6.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Pegasus
qpu_sampler = DWaveSampler(solver='Advantage2_prototype2.6',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Zephyr

# QUBO example
Q = {
    ('s1', 's1'): 1.0,
    ('s1', 's2'): 2.0,
    ('s1', 's3'): -2.0,
    ('s1', 's4'): -4.0,
    ('s2', 's2'): 1.0,
    ('s2', 's3'): -2.0,
    ('s2', 's4'): -4.0,
    ('s3', 's3'): 10.0,
    ('s3', 's4'): 4.0,
    ('s4', 's4'): 4.0,
}

# Reverse annealing schedule
ANNEAL_TIME = 20
anneal_schedule = [[0, 1], [ANNEAL_TIME / 2, 0.5], [ANNEAL_TIME, 1]
]

# Initial guess
initial_state ={'s1': 1, 's2': 1, 's3': 0, 's4': 1}
# Convert QUBO dictionary to a BinaryQuadraticModel
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

# Set up the sampler (using a D-Wave solver; adjust solver parameters as needed)
sampler = EmbeddingComposite(qpu_sampler)

# Solve the problem
sampleset = sampler.sample(
    bqm,
    num_reads=1000,
    initial_state=initial_state,
    anneal_schedule=anneal_schedule,
    reinitialize_state=False
)

print("Reverse Annealing Samples:")
for sample, energy in sampleset.data(['sample', 'energy']):
    print(sample, "Energy:", energy)


# sampleset.data(['sample', 'energy', 'num_occurrences']) can give us samples and their counts
# Convert samples to a tuple or string representation
solutions = []
for datum in sampleset.data(['sample', 'energy', 'num_occurrences']):
    sample = datum.sample
    # Convert dictionary (e.g. {'s1': 0, 's2': 0, 's3': 0}) to a string "000" or tuple (0,0,0)
    solution_str = ''.join(str(sample[var]) for var in sorted(sample.keys()))
    # Repeat each solution as many times as num_occurrences if you want to consider total frequency
    solutions.extend([solution_str] * datum.num_occurrences)

# Count how many times each solution appears
solution_counts = Counter(solutions)

# Prepare data for plotting
labels = list(solution_counts.keys())
counts = list(solution_counts.values())

# Create the histogram/bar chart
plt.figure(figsize=(8, 6))
plt.bar(labels, counts, color='navy')

plt.xlabel('Solutions')
plt.ylabel('Frequency')
plt.title('Histogram of QUBO Solutions')

# Optionally, display the counts above each bar
for i, count in enumerate(counts):
    plt.text(i, count + 0.1, str(count), ha='center', va='bottom')

plt.tight_layout()
plt.show()

# dwave.inspector.show(sampleset)
