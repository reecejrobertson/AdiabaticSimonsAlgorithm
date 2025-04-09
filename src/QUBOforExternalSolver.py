import numpy as np

MIN_N = 2
MAX_N = 50

def generateQUBO(n, penalties=None, matrix=False):
    '''
    Generate the QUBO for a Simon's problem. The "secret string" is the
    bitstring of all 1s.
    PARAMS:
        n       (int):   The number of qubits in the oracle.
        penalty (array): The penalties added to output transitions.
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
    if type(penalties) == None:
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
    
for n in range(MIN_N, MAX_N+1):

    with open(f'../QUBOs/balanced/qubo-{n}.txt', 'w+') as file:
        penalties = [(-1)**i * 2 for i in range(n)]
        file.write(str(generateQUBO(n, penalties=penalties)))

    with open(f'../QUBOs/random/qubo-{n}.txt', 'w+') as file:
        penalties = np.random.choice([2, -2], size=n, replace=True)
        file.write(str(generateQUBO(n, penalties=penalties)))