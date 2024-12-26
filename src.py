#TODO: For the moderately sized working problem:
    # TODO: 3d plot with annealing times [10ms, 20ms, 50ms, 100ms, 200ms, 500ms, 1000ms]
    # TODO: Try this from all starting states for reverse annealing
        # TODO: Add pausing into schedule
    # TODO: Repeat with forward annealing & fast annealing
    # TODO: 8192 samples up to 200ms, 1024 samples for 500-1000ms
    # TODO: Vary s3 variable

#---------------------------------------------------------------------------
#                            IMPORTS
#---------------------------------------------------------------------------

import dimod
import matplotlib.pyplot as plt
from collections import Counter
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

#---------------------------------------------------------------------------
#                            CONSTANTS
#---------------------------------------------------------------------------
VERBOSE = False
SHOTS = 1000
DPI = 300
N = 21

# QPU in Juelich
# qpu_sampler = DWaveSampler(solver='Advantage_system5.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4',region='eu-central-1')  # Pegasus Germany
# qpu_sampler = DWaveSampler(solver='Advantage_system4.1',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Pegasus 
# qpu_sampler = DWaveSampler(solver='Advantage_system6.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Pegasus
QPU_SAMPLER = DWaveSampler(solver='Advantage2_prototype2.6',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Zephyr

# Set up the sampler (using a D-Wave solver; adjust solver parameters as needed)
SAMPLER = EmbeddingComposite(QPU_SAMPLER)

# Reverse annealing schedule
ANNEAL_TIME = 20
ANNEAL_SCHEDULE = [[0, 1], [ANNEAL_TIME / 2, 0.5], [ANNEAL_TIME, 1]]

#---------------------------------------------------------------------------
#                            FUNCTIONS
#---------------------------------------------------------------------------
def generateQUBO(n, outputCouplingStrength=10):

    def indexIter():
        i = 1
        while True:
            yield f's{i}'
            i += 1

    Q = dict()
    qubitIndex = indexIter()
    in2 = next(qubitIndex)
    for i in range(n-1):
        in1 = in2
        in2 = next(qubitIndex)
        out = next(qubitIndex)
        ancilla = next(qubitIndex)
        Q[(in1, in1)] = 1.0
        Q[(in1, in2)] = 2.0
        Q[(in1, out)] = -2.0
        Q[(in1, ancilla)] = -4.0
        Q[(in2, in2)] = 1.0
        Q[(in2, out)] = -2.0
        Q[(in2, ancilla)] = -4.0
        Q[(out, out)] = 1.0 * outputCouplingStrength
        Q[(out, ancilla)] = 4.0
        Q[(ancilla, ancilla)] = 4.0
    return Q

def decodeBitstrings(bitstringList):
    decodedStrings = []
    for bitstring in bitstringList:
        bits = list(bitstring)
        inputs = []
        outputs = []
        ancillas = []
        inputs.append(bits.pop(0))
        while len(bits) > 0:
            inputs.append(bits.pop(0))
            outputs.append(bits.pop(0))
            ancillas.append(bits.pop(0))
        decodedStrings.append(''.join(inputs) + ' ' + ''.join(outputs) + ' ' + ''.join(ancillas))
    return decodedStrings

#---------------------------------------------------------------------------
#                            EXPERIMENT
#---------------------------------------------------------------------------

# For all n:
for n in range(2, N):

    print(f'n = {n}')

    # Generate a QUBO of size n
    Q = generateQUBO(n)

    # Convert QUBO dictionary to a BinaryQuadraticModel
    bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

    # For both annealing schedules:
    for schedule in ['Forward', 'Reverse']:

        print(f'{schedule} Annealing')

        if schedule == 'Reverse':

            # Initial guess
            initialState = dict()
            for i in range(1, 3*n - 1):
                initialState[f's{i}'] = 0

            # Solve the problem
            sampleSet = SAMPLER.sample(
                bqm,
                num_reads=SHOTS,
                initial_state=initialState,
                anneal_schedule=ANNEAL_SCHEDULE,
                reinitialize_state=False
            )

        else:

            # Solve the problem
            sampleSet = SAMPLER.sample(
                bqm,
                num_reads=SHOTS
            )

        # Print the sample results
        if VERBOSE:
            print(f'{schedule} Annealing Samples:')
            for sample, energy in sampleSet.data(['sample', 'energy']):
                print(sample, 'Energy:', energy)

        # sampleset.data(['sample', 'energy', 'num_occurrences']) can give us samples and their counts
        # Convert samples to a tuple or string representation
        solutions = []
        for datum in sampleSet.data(['sample', 'energy', 'num_occurrences']):
            sample = datum.sample

            # Convert dictionary (e.g. {'s1': 0, 's2': 0, 's3': 0}) to a string '000' or tuple (0,0,0)
            solutionStr = ''.join(str(sample[var]) for var in sorted(sample.keys()))

            # Repeat each solution as many times as num_occurrences if you want to consider total frequency
            solutions.extend([solutionStr] * datum.num_occurrences)

        # Count how many times each solution appears
        solutionCounts = Counter(solutions)

        # Prepare data for plotting
        labels = decodeBitstrings(list(solutionCounts.keys()))
        counts = list(solutionCounts.values())

        # Create the histogram/bar chart
        plt.figure(figsize=(8, 6))
        plt.bar(labels, counts, color='navy')
        
        # Label and title the figure
        plt.xlabel('Solutions')
        plt.ylabel('Frequency')
        plt.title(fr'Histogram of {schedule} Annealing Solutions for $n={n}$ QUBO')

        # Optionally, display the counts above each bar
        for i, count in enumerate(counts):
            plt.text(i, count + 0.1, str(count), ha='center', va='bottom')

        # Save the figure as an image
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'figs/{n}{schedule}.png', dpi=DPI)
        plt.close()

# dwave.inspector.show(sampleset)
