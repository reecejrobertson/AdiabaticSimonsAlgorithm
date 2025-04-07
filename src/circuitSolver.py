#-------------------------------------------------------------------------------
#                                 IMPORTS
#-------------------------------------------------------------------------------

import time
from qiskit import QuantumCircuit
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

#-------------------------------------------------------------------------------
#                                 CONSTANTS
#-------------------------------------------------------------------------------

SHOTS = 512
MAX_ORACLE_SIZE = 51
 
#-------------------------------------------------------------------------------
#                                 FUNCTIONS
#-------------------------------------------------------------------------------

def generateCircuit(n):
    circuit = QuantumCircuit(2*n, n)
    for m in range(n):
        circuit.h(m)
    for m in range(n):
        circuit.cx(m, n+m)
    for m in range(n):
        circuit.cx(0, n+m)
    for m in range(n):
        circuit.h(m)
    for m in range(n):
        circuit.measure(m, m)
    return circuit

#-------------------------------------------------------------------------------
#                                 OBJECTS
#-------------------------------------------------------------------------------

service = QiskitRuntimeService()
backend = service.least_busy(operational=True, simulator=False)
pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
sampler = Sampler(backend)
sampler.options.dynamical_decoupling.enable = True
sampler.options.dynamical_decoupling.sequence_type = 'XX'
print(backend.name)
print('-------------------------------------------------------------------')
 
#-------------------------------------------------------------------------------
#                                 EXPERIMENT
#-------------------------------------------------------------------------------

for n in range(2, MAX_ORACLE_SIZE):
    startTime = time.time()
    print(f'Oracle {n}')
    circuit = generateCircuit(n)
    isa_circuit = pm.run(circuit)
    print('Submitted')
    job = sampler.run([isa_circuit], shots=SHOTS)
    print('Done')
    results = job.result()[0]
    counts = results.data.c.get_counts()
    with open(f'../data/circuits/{backend.name}-{n}.txt', 'w+') as file:
        file.write(f'{str(counts)}\n')
        file.write(f'{str(job.usage_estimation)}\n')
        file.write(f'{str(isa_circuit.count_ops())}')
    print(f'Elapsed time: {time.time() - startTime}')
    print('-------------------------------------------------------------------')