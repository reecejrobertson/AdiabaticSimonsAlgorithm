import numpy as np

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

    if not matrix:

        # Return the final QUBO.
        return Q

    else:

        M = np.zeros((3*(n-1)+1, 3*(n-1)+1))
        # print(M)
        for var1, var2 in Q:
            M[var1, var2] = Q[(var1, var2)]
            M[var2, var1] = Q[(var1, var2)]

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

    # # Initialize an empty list of decoded results.
    # decodedStrings = []

    # # For each measured bitstring:
    # for bitstring in bitstringList:

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
    for i in range(len(outputs)):
        if (inputs[i] + inputs[i+1]) % 2 != outputs[i]:
            return False
    return True

def computeExpectation(n):
    Q = generateQUBO(n, 1, True)
    num_qubits = len(Q)
    temp = set()
    for i in range(2**num_qubits):
        x = np.array(list(map(int, list(format(i, f'0{num_qubits}b')))))
        # print(x)
        inputs, outputs, ancillas = decodeBitstring(x)
        if isInputValid(inputs, outputs, ancillas):
            key = ''.join(map(str, inputs)) + ' ' + ''.join(map(str, outputs))
            if key not in temp:
                temp.add(key)
                print(key)
        # print(decoded)
        # print(isInputValid(decoded[0], decoded[1], decoded[2]))
        # expectation = x.T @ Q @ x
        # print(expectation)
        # print(decodeBitstring(x))

computeExpectation(3)