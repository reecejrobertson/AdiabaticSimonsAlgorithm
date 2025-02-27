#---------------------------------------------------------------------------
#                            IMPORTS
#---------------------------------------------------------------------------

# import csv
import neal
import dimod
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite


#---------------------------------------------------------------------------
#                            CONSTANTS
#---------------------------------------------------------------------------

SIM = False                 # If true, then use D-Wave simulator.
VERBOSE = False             # If true, then logs are more detailed.
SHOTS = 512                 # The number of shots for each experiment.
DPI = 300                   # The resolution of figures.
ITERATIONS = range(1, 2)    # The number of iterations for all experiments.

if SIM:
    NS = range(2, 21)       # Range for number of qubits for simulator.
else:
    NS = [4, 8, 16, 32, 40]      # Range for number of qubits for hardware.

    # Get the API token for D-Wave.
    with open('../APIs/dwave.txt') as file:
        token = file.readline()

    # Access the D-Wave devices.
    ZEPHYR = DWaveSampler(solver='Advantage2_prototype2.6',token=token)
    # PEGASUS4_1 = DWaveSampler(solver='Advantage_system4.1',token=token)
    # PEGASUS6_4 = DWaveSampler(solver='Advantage_system6.4',token=token)
    QPU_SAMPLERS = {
        'zephyr': EmbeddingComposite(ZEPHYR),
        # 'pegasus6_4': EmbeddingComposite(PEGASUS6_4),
        # 'pegasus4_1': EmbeddingComposite(PEGASUS4_1)
}

# # The penalty parameters for output transitions.
# PENALTIES = [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]

# The annealing times.
# ANNEALING_TIMES = [10, 20, 50, 100, 200, 500, 1000]

# The annealing schedules.
ANNEALING_SCHEDULES = ['forward']#, 'reverse']

# The percentage of time to be spent pausing (reverse annealing only).
PAUSING_PERCENTAGES = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]


#---------------------------------------------------------------------------
#                            FUNCTIONS
#---------------------------------------------------------------------------

def generateQUBO(n, penalties=None, matrix=False):
    '''
    Generate the QUBO for a Simon's problem. The "secret string" is the
    bitstring of all 1s.
    PARAMS:
        n       (int): The number of qubits in the oracle.
        penalty (int): The penalty added to output transitions. None by default.
    RETURNS:
        (dict) The QUBO for the Simon's problem.
    '''

    def indexIter():
        '''
        An iterator to track qubit indices.
        '''
        i = 0
        while True:
            yield i
            i += 1

    # No penalties added by default.
    if penalties == None:
        penalties = [0] * (n - 1)

    # Initialize needed variables.
    Q = dict()
    qubitIndex = indexIter()
    in2 = next(qubitIndex)
    Q[(in2, in2)] = 0.0

    # For each qubit:
    for i in range(n-1):

        # Update qubit indices.
        in1 = in2
        in2 = next(qubitIndex)
        out = next(qubitIndex)
        ancilla = next(qubitIndex)

        # Add a one qubit to the QUBO oracle.
        Q[(in1, in1)] += 1.0
        Q[(in1, in2)] = 2.0
        Q[(in1, out)] = -2.0
        Q[(in1, ancilla)] = -4.0
        Q[(in2, in2)] = 1.0
        Q[(in2, out)] = -2.0
        Q[(in2, ancilla)] = -4.0
        Q[(out, out)] = 1.0 + penalties[i]
        Q[(out, ancilla)] = 4.0
        Q[(ancilla, ancilla)] = 4.0

    if not matrix:

        # Return the final QUBO.
        return Q

    else:

        M = np.zeros((3*(n-1)+1, 3*(n-1)+1))
        for var1, var2 in Q:
            M[var1, var2] = Q[(var1, var2)]

        print(M)
        return M

def decodeBitstring(bitstringList):
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
    labels = decodeBitstring(list(resultCounts.keys()))
    counts = list(resultCounts.values())

    # Create a histogram/bar chart.
    plt.figure(figsize=(8, 6))
    plt.tight_layout()
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
    plt.savefig(f'../figs/specialPenalties/{title}.png', dpi=DPI)
    plt.close()

    # Write the data to a file.
    with open(f'../data/specialPenalties/{title}.txt', 'w+') as file:
        for i, label in enumerate(labels):
            file.write(f'{label}:{counts[i]}\n')

def generatePenalties(n, p=1, option=1):
    penalties = []
    if option == 1:
        for i in range(n-1):
            if i % 2 == 0:
                penalties.append(p)
            else:
                penalties.append(-p)
    if option == 2:
        penalties += [p] * int(n/2)
        penalties += [-p] * ((n-1) - int(n/2))
    if option == 3:
        mask = np.random.choice(list(range(n-1)), int(n/2))
        for i in range(n-1):
            if i in mask:
                penalties.append(-p)
            else:
                penalties.append(p)
    return penalties


#---------------------------------------------------------------------------
#                            EXPERIMENT
#---------------------------------------------------------------------------

# For each iteration:
for i in ITERATIONS:
    print(f'i = {i}')

    if SIM:
        print("simulating")

        # For all n:
        for n in NS:
            print(f'n = {n}')

            penalties = [[0]*(n-1), [1]*(n-1), [n]*(n-1)]

            # For each penalty parameter:
            for penalty in penalties:
                print(f'penalty = {penalty}')

                # Generate a QUBO of size n.
                Q = generateQUBO(n, penalty)

                # Convert QUBO dictionary to a BinaryQuadraticModel.
                bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

                # Get the simulated sampler.
                s = neal.SimulatedAnnealingSampler()

                # Set sweeping parameters by intuition.
                num_sweeps = 1000
                beta_schedule_type="geometric"
                beta_range = (0.01, 100)
        
                # Generate the title for the experiment.
                title = (
                    "SimulatedAnnealingSampler" + '-' +
                    str(n)              + '-' +
                    str(penalty[0])  + '-' +
                    str(num_sweeps)     + '-' +
                    str(beta_schedule_type)
                )

                # Run the experiment.
                sampleSet = s.sample(
                    bqm,
                    num_reads=SHOTS,
                    beta_range = beta_range, 
                    num_sweeps = num_sweeps,
                    beta_schedule_type=beta_schedule_type
                    )

                # Extract and format the experimental resuls.
                extractResults(sampleSet, title)

    else:

        # For each D-Wave device:
        for qpuSamplerName in QPU_SAMPLERS.keys():
            print(f'device = {qpuSamplerName}')
            qpuSampler = QPU_SAMPLERS[qpuSamplerName]

            # For all n:
            for n in NS:
                print(f'n = {n}')

                # penalties = [[2]*(n-1), [5]*(n-1), [n/2]*(n-1)]
                penalties = [
                    generatePenalties(n, 2, 1),
                    generatePenalties(n, 2, 2),
                    generatePenalties(n, 2, 3),
                    generatePenalties(n, 5, 1),
                    generatePenalties(n, 5, 2),
                    generatePenalties(n, 5, 3),
                ]
                annealingTimes = [100, 1000, n*10]

                # For each penalty parameter:
                for penalty in penalties:
                    print(f'penalty = {penalty}')

                    # Generate a QUBO of size n.
                    Q = generateQUBO(n, penalty)

                    # Convert QUBO dictionary to a BinaryQuadraticModel.
                    bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

                    # For each annealing time:
                    for time in annealingTimes:
                        print(f'time = {time}')

                        # For both annealing schedules:
                        for schedule in ANNEALING_SCHEDULES:
                            print(f'annealing schedule = {schedule}')

                            # If the schedule is the reserse annealing:
                            if schedule == 'reverse':

                                # Then for each pausing percentage:
                                for pausingPercentage in PAUSING_PERCENTAGES:
                                    print(
                                        f'pausing % = {pausingPercentage}'
                                    )

                                    # Generate the title of the experiment.
                                    title = (
                                        str(qpuSamplerName) + '-' +
                                        str(n)              + '-' +
                                        str(penalty[0])     + '-' +
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

                                    # Extract and format the results.
                                    extractResults(sampleSet, title)

                            # If the scheulde is forward annealing:
                            else:

                                # Generate the title for the experiment.
                                title = (
                                    str(qpuSamplerName) + '-' +
                                    str(n)              + '-' +
                                    str(penalty)     + '-' +
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
