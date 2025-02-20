#---------------------------------------------------------------------------
#                            IMPORTS
#---------------------------------------------------------------------------

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

VERBOSE = False             # If true, logs are more detailed.
SHOTS = 512                 # The number of shots for each experiment.
DPI = 300                   # The resolution of figures.
ITERATIONS = range(1, 2)    # The number of iterations for all experiments.
NS = [3]                    # Range for number of qubits for experiments.

# The penalty parameters for output transitions.
# PENALTIES = [0.1, 1, 2, 5, 10, 20, 50, 100]
# PENALTIES = [-10**2, -10**3, -10**4, -10**5, -10**6, -10**7, -10**8]
PENALTIES = [0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]

#---------------------------------------------------------------------------
#                            FUNCTIONS
#---------------------------------------------------------------------------

def generateQUBO(n, penalties=None, matrix=False):
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
        i = 0
        while True:
            yield i
            i += 1

    # Set the default penalty parameter to scale proportionally to oracle size.
    if penalties == None:
        penalties = [2**n] * (n - 1)

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
        Q[(out, out)] = 1.0 * penalties[i]
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
    plt.savefig(f'simFigs/{title}.png', dpi=DPI)
    plt.close()

    # Write the data to a file.
    with open(f'simData/{title}.txt', 'w+') as file:
        for i, label in enumerate(labels):
            file.write(f'{label}:{counts[i]}\n')

#---------------------------------------------------------------------------
#                            EXPERIMENT
#---------------------------------------------------------------------------

# For each iteration:
for i in ITERATIONS:
    print(f'i = {i}')
    print("simulating")

    # For all n:
    for n in NS:
        print(f'n = {n}')

        # For each penalty parameter:
        for penalty1 in PENALTIES:
            print(f'first penalty =  {penalty1}')

            for penalty2 in PENALTIES:
                print(f'second penalty = {penalty2}')

                penalties = [penalty1, penalty2]

                # Generate a QUBO of size n.
                Q = generateQUBO(n, penalties)

                # Convert QUBO dictionary to a BinaryQuadraticModel.
                bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

                # Run the experiment.
                s = neal.SimulatedAnnealingSampler()
                # https://docs.ocean.dwavesys.com/projects/neal/en/latest/reference/generated/neal.sampler.SimulatedAnnealingSampler.sample.html

                # these 3 are sweeping parameters - I set them by intuition
                num_sweeps = 1000
                beta_schedule_type="geometric"
                beta_range = (0.01, 100)
        
                # Generate the title for the experiment.
                title = (
                    "SimulatedAnnealingSampler" + '-' +
                    str(n)              + '-' +
                    str('_'.join(map(str, penalties)))        + '-' +
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
