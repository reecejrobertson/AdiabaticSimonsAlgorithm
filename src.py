import dimod
import matplotlib.pyplot as plt
from collections import Counter
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

# QPU in Juelich
qpu_sampler = DWaveSampler(solver='Advantage_system5.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4',region="eu-central-1")

# QUBO example
Q = {
    ('s1', 's1'): 1.0,
    ('s2', 's2'): 1.0,
    ('s3', 's3'): 2.0,
    ('s1', 's3'): -2.0,
    ('s2', 's3'): -2.0
}

# Reverse annealing schedule
ANNEAL_TIME = 200
anneal_schedule = [[0, 1], [ANNEAL_TIME / 2, 0.5], [ANNEAL_TIME, 1]
]

# Initial guess
initial_state ={'s1': 1, 's2': 1, 's3': 0}
# Convert QUBO dictionary to a BinaryQuadraticModel
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

# Set up the sampler (using a D-Wave solver; adjust solver parameters as needed)
sampler = EmbeddingComposite(qpu_sampler)

# Solve the problem
sampleset = sampler.sample(
    bqm,
    num_reads=100,
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
