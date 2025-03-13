# Adiabatic Simon's Algorithm

This repository contains an implementation of an **adiabatic version of Simon's Algorithm**.
The algorithm is a quantum algorithm that solves the **Simon's problem**, which can be solved exponentially faster using quantum resources than classical ones.

## Table of Contents

- [Introduction](#introduction)
- [Algorithm](#algorithm)
- [Structure](#structure)
- [Usage](#usage)
- [License](#license)

## Introduction

Simon's problem involves finding a hidden string $\( s \)$ such that a black-box function $\( f(x) \)$ satisfies:
- $\( f(x) = f(x') \)$ if $\( x \oplus x' = s \)$
- $\( f(x) \neq f(x') \)$ otherwise.

This implementation provides an adiabatic version of Simon's algorithm which uses quantum adiabatic evolution to find the hidden string $\( s \).$
For more information on the theory and implementation details, see the research paper—link forthcoming.

## Algorithm

The adiabatic version of Simon's algorithm adapts the problem to a quantum adiabatic framework, which minimizes the ground state energy of a Hamiltonian whose ground state encodes the solution to the problem.
In practice, we accomplish this by encoding an instance of Simon's problem in a QUBO which includes input, output, and ancillary qubits.
We enforce a penalty parameter on transitions of the output qubits, which effectively fixes the output.
We then allow the system to evolve, and recover the input pair that matches the fixed output.

## Structure

The repository structure is as follows:

```
adiabatic-simons-algorithm/
├── src/                    # Source code for the adiabatic Simon's algorithm
│   ├── analystcalSolver.py # Solve Simon's QUBO exactly
│   ├── dataProcessing.py   # Process generated data
│   ├── dwaveSolver.py      # Solve Simon's QUBO on D-Wave hardware
│   └── simulatorSolver.py  # Solve Simon's QUBO on D-Wave simulators
├── APIs/                   # API folder (must be created locally)
│   ├── dwave.py            # D-Wave API file (must be created locally)
├── data/                   # Data folder
├── archive/                # Old code files
├── .gitignore              # Git ignore file (e.g., system files)
├── LICENSE                 # License file
└── README.md               # Project description and instructions
```

Of these files, those in the source directory are the most important.
In particular, the `generateQUBO()` function in the `src/dwaveSolver.py` file generates QUBOs for Simon's problem of varying sizes, with varying penalty parameters.
Early experimental results indicate that a size $n$ penalty array with alternating entries $+2$ and $-2$ works well in most instances.

## Usage

To run the experiments on a D-Wave machine, provide your API token in the `APIs/dwave.txt` file, then run `python src/dwaveSolver.py`.
To run the experiments on a D-Wave simulator, run `python src/simulatorSolver.py`.
To process data obtained through experiments (including the data already contained in the `data/` directory), run `python src/dataProcessing.py`.
To compute the analytical solution to the problem, run `python src/analyticalSolver.py`

## License

This project is licensed under the MIT License - see the LICENSE file for details.
