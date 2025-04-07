import csv
import numpy as np
from matplotlib import pyplot as plt

def generateQUBO(n, penalty=None, matrix=False):
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
    if penalty == None:
        penalty = -2**n

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
        Q[(out, out)] = 1.0 * (penalty**(i+1))
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
    
def decodeBitstring(bitstring):
    '''
    Decode the results of many shots of an instance of Simon's problem.
    This involves identifying which qubits correspond to the input, output, and
    ancilla registers
    PARAMS:
        bitstringsList (list): A list of the measured bitstrings.
    RETURNS:
        (list): The decoded results.
    '''

    # Initialize an array for the input, output, and ancilla registers.
    inputs = []
    outputs = []
    ancillas = []

    # Extract the bits of the bitstring that relate to each register.
    bits = bitstring.tolist()
    inputs.append(bits.pop(0))
    while len(bits) > 0:
        inputs.append(bits.pop(0))
        outputs.append(bits.pop(0))
        ancillas.append(bits.pop(0))
    
    # Return the decoded results.
    return inputs, outputs, ancillas

def isInputValid(inputs, outputs, ancillas):
    '''
    Determines if a given input, output, ancilla tuple is a valid oracle output.
    PARAMS:
        inputs (array):   The bits of the input register.
        outputs (array):  The bits of the output register.
        ancillas (array): The bits of the ancilla register.
    RETURNS:
        (bool): True of the bits tuple is valid, otherwise False.
    '''

    # For each bit:
    for i in range(len(outputs)):

        # Check that the output and oracle constraints are satisfied.
        if (
            ((inputs[i] + inputs[i+1]) % 2 != outputs[i]) or
            ((inputs[i] and inputs[i+1]) != ancillas[i])
        ):
            return False
    
    # Return true only if all constraints are satisfied.
    return True

def computeExpectation():
    '''
    Analytically compute the expectation of the function of size 3.
    '''

    # Set constant variables.
    N = 3
    PENALTIES = [0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]
    
    # Instantiate a figure and a csv file.
    fig, axs = plt.subplots(3, 4, tight_layout=True, figsize=(12, 9))
    csvHeader = ['bitstring', 'expectation', 'penalty', 'valid']
    with open(f'data{N}.csv', 'w+', encoding='UTF8') as file:
        writer = csv.writer(file)
        writer.writerow(csvHeader)
    
    # For each penalty value:
    for idx, penalty in enumerate(PENALTIES):

        # Generate the corresponding QUBO.
        Q = generateQUBO(N, penalty, True)

        # Get the number of qubits in the problem.
        num_qubits = len(Q)

        # Instantiate needed data structures.
        valid = dict()
        invalid = []
        csvData = []

        # For each possible assignment of classical values to the qubits:
        for i in range(2**num_qubits):

            # Compute the expectation of that bitstring with the QUBO.
            x = np.array(list(map(int, list(format(i, f'0{num_qubits}b')))))
            expectation = x.T @ Q @ x

            # Decode the bitstring and check if it is valid.
            inputs, outputs, ancillas = decodeBitstring(x)
            isValid = isInputValid(inputs, outputs, ancillas)

            # Record the bitstring and its expectation in the proper structures.
            if isValid:
                key = ''.join(map(str, outputs))
                if key in valid:
                    valid[key].append(expectation)
                else:
                    valid[key] = [expectation]
            else:
                invalid.append(expectation)
            decodedString = str(
                ''.join(map(str, inputs)) + ' ' +
                ''.join(map(str, outputs)) + ' ' +
                ''.join(map(str, ancillas))
            )
            csvLine = [decodedString, expectation, penalty, isValid]
            csvData.append(csvLine)

        # Plot the results in a histogram.
        categories = [valid[result] for result in valid]
        categories.append(invalid)
        axs[int(idx/4), int(idx%4)].hist(
            categories,
            bins=100,
            color=['r', 'b', 'g', 'm', 'k'],
            histtype='barstacked',
            label=['Valid', 'Valid', 'Valid', 'Valid', 'Invalid']
        )

        # Label the plot.
        axs[int(idx/4), int(idx%4)].set_title(f'$p={penalty}$')
        if idx == 3:
            axs[int(idx/4), int(idx%4)].legend(loc='best')
        if idx > 7:
            axs[int(idx/4), int(idx%4)].set_xlabel(f'Energy')
        if idx in [0, 4, 8]:
            axs[int(idx/4), int(idx%4)].set_ylabel(f'Count')

        # Write the csv file.
        with open(f'data{N}.csv', 'a', encoding='UTF8') as file:
            writer = csv.writer(file)
            writer.writerows(csvData)

    # Save the plot.
    plt.savefig(f'fig{N}.png', dpi=300)

# Compute the expectation for the 3 qubit oracle.
computeExpectation()