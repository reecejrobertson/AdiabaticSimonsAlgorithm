# Adiabatic Simon's Algorithm

This repository contains an implementation of an **adiabatic version of Simon's Algorithm**.
A complete description of the algorithm can be found in the paper associated with this project: https://arxiv.org/abs/2504.10771.

## Table of Contents

- [Introduction](#introduction)
- [Algorithm](#algorithm)
- [Structure](#structure)
- [Usage](#usage)
- [License](#license)

## Introduction

Let $\\{0,1\\}^n$ be the set of all binary bitstrings of length $n$.
Suppose that one is given a black-box oracle function $f(x):\\{0,1\\}^n\rightarrow\\{0,1\\}^{n-1}$ with the property that there exists some fixed $s\in\\{0,1\\}^n$ such that $f(x) = f(x')$ for all $x \neq x'$ if and only if $x \oplus x' = s$.
In other words, $f(x)$ is a two-to-one periodic function with period $s$.
Simon's problem is to identify $s$, where the only allowable operation is to query $f$.
Simon showed that a fault-tolerant gate-based quantum computer can solve this problem with exponential advantage over a classical computer (https://epubs.siam.org/doi/10.1137/S0097539796298637)

This repository provides an adiabatic implementation of Simon's algorithm, which uses quantum adiabatic evolution to find the hidden string $s$.
At a high level, the algorithm works by encoding an instance of Simon's problem into a Quadratic Unconstrained Binary Optimization (QUBO) problem, which can be solved on a quantum annealer.
Once a QUBO has been obtained, penalty parameters are introduced on the output qubits.
Doing this prepares a degenerate ground state which contains two valid evaluations of the oracle which share an output, that is, our ground state contains both $z$ and $z'=z\oplus s$ for some fixed $z\in\\{0,1\\}^n$.
Assuming that both $z$ and $z'$ are sampled from the ground state with equal probability, then successfully sampling this state a few times is sufficient to solve the problem with high probability.
For the full details of the implementation and performance of the algorithm, refer to https://arxiv.org/abs/2504.10771.

## Structure

The repository structure is as follows:

```
adiabatic-simons-algorithm/
├── src/                              # Source code for the adiabatic Simon's algorithm
│   ├── analyticalSolver.py           # Solve Simon's QUBO exactly
│   ├── circuitSolver.py              # Solve Simon's problem via gate-based circuit implementation
│   ├── dwaveEmbedding.py             # View the embedding of a QUBO on a D-Wave device
│   ├── dwaveProcessing.py            # Process data generated through dwaveSolver.py
│   ├── dwaveSolver.py                # Solve Simon's QUBO on D-Wave hardware or simulator
│   ├── formatDataForPlotting.ipynb   # Format data for plotting with generatePlots.ipynb
│   ├── generatePlots.ipynb           # Generate clean plots to save in the figs/ folder
│   └── QUBOforExternalSolver.py      # Generate QUBOs to offload to external solvers
├── APIs/                             # API folder (must be created locally)
│   ├── dwave.py                      # D-Wave API file (must be created locally)
├── data/                             # Data folder
├── figs/                             # Figure folder
├── QUBOs/                            # Folder containing text files for balanced and random QUBOS
├── archive/                          # Old code, data, and figure files
├── .gitignore                        # Git ignore file (e.g., system files)
├── LICENSE                           # License file
└── README.md                         # Project description and instructions
```

Of these files, those in the source directory are the most important.
In particular, the `generateQUBO()` function in the `src/QUBOforExternalSolver.py` file generates QUBOs for Simon's problem of varying sizes and penalty parameters.
These QUBOs are saved in the `QUBOs/` directory and can be imported into external solvers.
Experimental results indicate that a size $n$ penalty array with alternating entries $+2$ and $-2$ works well in most instances.

## Usage

To run the experiments on a D-Wave machine, provide your API token in the `APIs/dwave.txt` file, then run `python src/dwaveSolver.py`.
To run the experiments on a D-Wave simulator, run `python src/simulatorSolver.py`.
To process data obtained through experiments (including the data already in the `data/` directory), run `python src/dataProcessing.py`.
To compute the analytical solution to the problem, run `python src/analyticalSolver.py`

## License

This project is licensed under the MIT License - see the LICENSE file for details.
