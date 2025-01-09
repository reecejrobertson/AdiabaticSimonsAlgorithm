#TODO: For the moderately sized working problem:
    # DONE: 3d plot with annealing times [10ms, 20ms, 50ms, 100ms, 200ms, 500ms, 1000ms]
    # DONE: Try this from all starting states for reverse annealing
        # DONE: Add pausing into schedule
    # TODO: Repeat with forward annealing & fast annealing
    # DONE: 8192 samples up to 200ms, 1024 samples for 500-1000ms
    # DONE: Vary s3 variable
    # TODO: Experiment with simulated annealing
    # TODO: Test QUBO on a heuristic simulator
    # TODO: Penalty inversely proportional to the system size

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

# Define the needed constants
VERBOSE = False
SHOTS = 512
DPI = 300

# Define the number of iterations
ITERATIONS = range(1, 2)

# Set up the sampler (using a D-Wave solver; adjust solver parameters as needed)
PEGASUS5_4 = DWaveSampler(solver='Advantage_system5.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4',region='eu-central-1')  # Pegasus Germany
PEGASUS4_1 = DWaveSampler(solver='Advantage_system4.1',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Pegasus 
PEGASUS6_4 = DWaveSampler(solver='Advantage_system6.4',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Pegasus
ZEPHYER = DWaveSampler(solver='Advantage2_prototype2.6',token='julr-a86ece088ec3ae431ae7ee0541c03112c43d7af4')  # Zephyr
QPU_SAMPLERS = {
    # 'pegasus5_4': EmbeddingComposite(PEGASUS5_4),
    'zephyer': EmbeddingComposite(ZEPHYER),
    'pegasus6_4': EmbeddingComposite(PEGASUS6_4),
    'pegasus4_1': EmbeddingComposite(PEGASUS4_1)
}

# Define the problem sizes
NS = [7,8]

# Define the output coupling strengths
OUTPUT_COUPLING_STRENGTHS = [1, 2, 5, 10, 20, 50, 100]

# Define the annealing times
ANNEALING_TIMES = [10, 20, 50, 100, 200, 500, 1000]

# Define the annealing schedules
ANNEALING_SCHEDULES = ['forward', 'reverse']

# Define pausing parameters for the reverse annealing schedule
PAUSING_PERCENTAGES = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

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

def generateReverseAnnealingSchedule(annealingTime, pausingPercentage):
    if pausingPercentage != 0:
        activePercentage = 1 - pausingPercentage
        pauseStart = activePercentage / 2. * annealingTime
        pauseEnd = pauseStart + (annealingTime * pausingPercentage)
        annealingSchedule = [[0, 1], [pauseStart, 0.5], [pauseEnd, 0.5], [annealingTime, 1]]
    else:
        annealingSchedule = [[0, 1], [annealingTime / 2., 0.5], [annealingTime, 1]]
    return annealingSchedule

def runExperiment(sampleSet, title):
    # Print the sample results
    if VERBOSE:
        print(f'annealing samples:')
        for sample, energy in sampleSet.data(['sample', 'energy']):
            print(sample, 'energy:', energy)

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

    # Optionally, display the counts above each bar
    for i, count in enumerate(counts):
        plt.text(i, count + 0.1, str(count), ha='center', va='bottom')

    # Save the figure as an image
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'figs/{title}.png', dpi=DPI)
    plt.close()

    with open(f'data/{title}.txt', 'w+') as file:
        for i, label in enumerate(labels):
            file.write(f'{label}:{counts[i]}\n')

#---------------------------------------------------------------------------
#                            EXPERIMENT
#---------------------------------------------------------------------------

for i in ITERATIONS:
    print(f'i = {i}')

    for qpuSamplerName in QPU_SAMPLERS.keys():
        print(f'device = {qpuSamplerName}')

        qpuSampler = QPU_SAMPLERS[qpuSamplerName]

        # For all n:
        for n in NS:
            print(f'n = {n}')

            for outputCouplingStrength in OUTPUT_COUPLING_STRENGTHS:
                print(f'output coupling strength = {outputCouplingStrength}')

                # Generate a QUBO of size n
                Q = generateQUBO(n, outputCouplingStrength)

                # Convert QUBO dictionary to a BinaryQuadraticModel
                bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

                for time in ANNEALING_TIMES:
                    print(f'time = {time}')

                    # For both annealing schedules:
                    for schedule in ANNEALING_SCHEDULES:
                        print(f'annealing schedule = {schedule}')

                        if schedule == 'reverse':

                            for pausingPercentage in PAUSING_PERCENTAGES:
                                print(f'pausing percentage = {pausingPercentage}')

                                title = f'{qpuSamplerName}-{n}-{outputCouplingStrength}-{time}-{schedule}-{pausingPercentage}'

                                annealingSchedule = generateReverseAnnealingSchedule(time, pausingPercentage)

                                # Initial guess
                                initialState = dict()
                                for i in range(1, 3*n - 1):
                                    initialState[f's{i}'] = 0

                                # Solve the problem
                                sampleSet = qpuSampler.sample(
                                    bqm,
                                    num_reads=SHOTS,
                                    initial_state=initialState,
                                    anneal_schedule=annealingSchedule,
                                    reinitialize_state=False
                                )

                                runExperiment(sampleSet, title)

                        else:

                            title = f'{qpuSamplerName}-{n}-{outputCouplingStrength}-{time}-{schedule}'

                            # Solve the problem
                            sampleSet = qpuSampler.sample(
                                bqm,
                                num_reads=SHOTS,
                                annealing_time = time
                            )

                            runExperiment(sampleSet, title)
