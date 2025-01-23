#---------------------------------------------------------------------------
#                             TODO LIST
#---------------------------------------------------------------------------

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

# TODO: Generalize this to larger problem [goal: n=20]

#---------------------------------------------------------------------------
#                             SIMULATOR CODE
#---------------------------------------------------------------------------

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

VERBOSE = False             # If true, logs are more detailed.
SHOTS = 512                 # The number of shots for each experiment.
DPI = 300                   # The resolution of figures.
ITERATIONS = range(1, 2)    # The number of iterations for all experiments.
NS = [3, 8]                 # Range for number of qubits for experiments.

# The penalty parameters for output transitions.
PENALTIES = [1, 2, 5, 10, 20, 50, 100]

# The annealing times.
ANNEALING_TIMES = [10, 20, 50, 100, 200, 500, 1000]

# The annealing schedules.
ANNEALING_SCHEDULES = ['forward', 'reverse']

# The percentage of time to be spent pausing (reverse annealing only).
PAUSING_PERCENTAGES = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

# Get the API token for D-Wave.
with open('APIs/dwave.txt') as file:
    token = file.readline()

# Access the D-Wave devices.
PEGASUS4_1 = DWaveSampler(solver='Advantage_system4.1',token=token) 
PEGASUS6_4 = DWaveSampler(solver='Advantage_system6.4',token=token)
ZEPHYER = DWaveSampler(solver='Advantage2_prototype2.6',token=token)
QPU_SAMPLERS = {
    'zephyer': EmbeddingComposite(ZEPHYER),
    'pegasus6_4': EmbeddingComposite(PEGASUS6_4),
    'pegasus4_1': EmbeddingComposite(PEGASUS4_1)
}

#---------------------------------------------------------------------------
#                            FUNCTIONS
#---------------------------------------------------------------------------

def generateQUBO(n, penalty=10):
    '''
    Generate the QUBO for a Simon's problem. The "secret string" is the
    bitstring of all 1s.
    PARAMS:
        n       (int):  The number of qubits in the oracle.
        penalty (int):  The penalty applied to output transitions.
    RETURNS:
        (dict) The QUBO for the Simon's problem.
    '''

    def indexIter():
        '''
        An iterator to track qubit indices.
        '''
        i = 1
        while True:
            yield f's{i}'
            i += 1

    # Initialize needed variables.
    Q = dict()
    qubitIndex = indexIter()
    in2 = next(qubitIndex)

    # For each qubit:
    for i in range(n-1):

        # Update qubit indices.
        in1 = in2
        in2 = next(qubitIndex)
        out = next(qubitIndex)
        ancilla = next(qubitIndex)

        # Add a one qubit to the QUBO oracle.
        Q[(in1, in1)] = 1.0
        Q[(in1, in2)] = 2.0
        Q[(in1, out)] = -2.0
        Q[(in1, ancilla)] = -4.0
        Q[(in2, in2)] = 1.0
        Q[(in2, out)] = -2.0
        Q[(in2, ancilla)] = -4.0
        Q[(out, out)] = 1.0 * penalty
        Q[(out, ancilla)] = 4.0
        Q[(ancilla, ancilla)] = 4.0

    # Return the final QUBO.
    return Q

def decodeBitstrings(bitstringList):
    '''
    Decode the results of many shots of an instance of Simon's problem.
    This involves identifying which qubits correspond to the input, output, and
    ancilla registers
    PARAMS:
        bitstringsList (list): A list of the measured bitstrings.
    RETURNS:
        (list): The decoded results.
    '''

    # Initialize an empty list of decoded results.
    decodedStrings = []

    # For each measured bitstring:
    for bitstring in bitstringList:

        # Initialize an array for the input, output, and ancilla registers.
        inputs = []
        outputs = []
        ancillas = []

        # Extract the bits of the bitstring that relate to each register.
        bits = list(bitstring)
        inputs.append(bits.pop(0))
        while len(bits) > 0:
            inputs.append(bits.pop(0))
            outputs.append(bits.pop(0))
            ancillas.append(bits.pop(0))
        
        # Append the decoded result to the list of decoded results.
        decodedStrings.append(
            ''.join(inputs)  + ' ' +
            ''.join(outputs) + ' ' +
            ''.join(ancillas)
        )
    
    # Return the decoded results.
    return decodedStrings

def setReverseSchedule(annealingTime, pausingPercentage):
    '''
    Generates a schedule for a reverse annealing experiment.
    PARAMS:
        annealingTime       (float): The time to run the annealer.
        pausingPercentage   (float): The percentage of the runtime to pause.
    RETURNS:
        (list): The schedule for the annealer.
    '''

    # If there is a nonzero pausing percentage:
    if pausingPercentage != 0:

        # Compute the time to begin and end the pause.
        activePercentage = 1 - pausingPercentage
        pauseStart = activePercentage / 2. * annealingTime
        pauseEnd = pauseStart + (annealingTime * pausingPercentage)
    
        # Define the schedule.
        annealingSchedule = [
            [0, 1],
            [pauseStart, 0.5],
            [pauseEnd, 0.5],
            [annealingTime, 1]
        ]
    
    # Otherwise:
    else:

        # Define the schedule without a pause.
        annealingSchedule = [
            [0, 1],
            [annealingTime / 2., 0.5],
            [annealingTime, 1]
        ]
    
    # In both cases, return the schedule.
    return annealingSchedule

def extractResults(jobResults, title):
    '''
    Run an instance of the experiment with specified variables. Plots and saves
    a figure representing the data results, and saves the raw data to a file.
    PARAMS:
        jobResults  (Future): The results of the D-Wave job.
        title       (string): The title for the figure.
    '''
    
    # Initialize an empty list for the measurement results.
    results = []

    # If versose, print the sample results.
    if VERBOSE:
        print(f'annealing samples:')
        for sample, energy in jobResults.data(['sample', 'energy']):
            print(sample, 'energy:', energy)

    # For each measurement result:
    for datum in jobResults.data(['sample', 'energy', 'num_occurrences']):

        # Extract the measurement results and convert to a string.
        sample = datum.sample
        solutionStr = ''.join(str(sample[var]) for var in sorted(sample.keys()))

        # Repeat each result according to its number of occurances.
        results.extend([solutionStr] * datum.num_occurrences)

    # Count the number of instances of each result.
    resultCounts = Counter(results)

    # Prepare data for plotting.
    labels = decodeBitstrings(list(resultCounts.keys()))
    counts = list(resultCounts.values())

    # Create a histogram/bar chart.
    plt.figure(figsize=(8, 6))
    plt.bar(labels, counts, color='navy')
    
    # Label and title the figure.
    plt.xlabel('Solutions')
    plt.ylabel('Frequency')

    # Display the counts above each bar (optional).
    for i, count in enumerate(counts):
        plt.text(i, count + 0.1, str(count), ha='center', va='bottom')

    # Save the figure as an image.
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'figs/{title}.png', dpi=DPI)
    plt.close()

    # Write the data to a file.
    with open(f'data/{title}.txt', 'w+') as file:
        for i, label in enumerate(labels):
            file.write(f'{label}:{counts[i]}\n')

#---------------------------------------------------------------------------
#                            EXPERIMENT
#---------------------------------------------------------------------------

# For each iteration:
for i in ITERATIONS:
    print(f'i = {i}')

    # For each D-Wave device:
    for qpuSamplerName in QPU_SAMPLERS.keys():
        print(f'device = {qpuSamplerName}')
        qpuSampler = QPU_SAMPLERS[qpuSamplerName]

        # For all n:
        for n in NS:
            print(f'n = {n}')

            # For each penalty parameter:
            for penalty in PENALTIES:
                print(f'penalty = {penalty}')

                # Generate a QUBO of size n.
                Q = generateQUBO(n, penalty)

                # Convert QUBO dictionary to a BinaryQuadraticModel.
                bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

                # For each annealing time:
                for time in ANNEALING_TIMES:
                    print(f'time = {time}')

                    # For both annealing schedules:
                    for schedule in ANNEALING_SCHEDULES:
                        print(f'annealing schedule = {schedule}')

                        # If the schedule is the reserse annealing:
                        if schedule == 'reverse':

                            # Then for each pausing percentage:
                            for pausingPercentage in PAUSING_PERCENTAGES:
                                print(
                                    f'pausing percentage = {pausingPercentage}'
                                )

                                # Generate the title of the experiment.
                                title = (
                                    str(qpuSamplerName) + '-' +
                                    str(n)              + '-' +
                                    str(penalty)        + '-' +
                                    str(time)           + '-' + 
                                    str(schedule)       + '-' +
                                    str(pausingPercentage)
                                )

                                # Deterime the annealing schedule.
                                annealingSchedule = setReverseSchedule(
                                    time,
                                    pausingPercentage
                                )

                                # Generate an initial guess.
                                initialState = dict()
                                for i in range(1, 3*n - 1):
                                    initialState[f's{i}'] = 0

                                # Run the experiment.
                                sampleSet = qpuSampler.sample(
                                    bqm,
                                    num_reads=SHOTS,
                                    initial_state=initialState,
                                    anneal_schedule=annealingSchedule,
                                    reinitialize_state=False
                                )

                                # Extract and format the experimental results.
                                extractResults(sampleSet, title)

                        # If the scheulde is forward annealing:
                        else:

                            # Generate the title for the experiment.
                            title = (
                                str(qpuSamplerName) + '-' +
                                str(n)              + '-' +
                                str(penalty)        + '-' +
                                str(time)           + '-' +
                                str(schedule)
                            )

                            # Run the experiment.
                            sampleSet = qpuSampler.sample(
                                bqm,
                                num_reads=SHOTS,
                                annealing_time = time
                            )

                            # Extract and format the experimental resuls.
                            extractResults(sampleSet, title)
