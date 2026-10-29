# Diagonalization of symmetric tridiagonal matrix (quantum harmonic oscillator problem).

The purpouse of the project is to solve a 1D quantum harmonic oscillator problem numerically. The second derivative is replaced with a three-point finite difference formula on a uniform grid x ∈ [-x_max, x_max],
which gives a symmetric tridiagonal Hamiltonian matrix. The matrix is
diagonalized with LAPACKE (`LAPACKE_dstev`). Numerical eigenvalues and
eigenvectors are compared with the exact solutions.

**Author:** Michal Sikora
**Subject:** Numerical Methods in Physics
**Major and Faculty:** Computer Physics, Faculty of Physics and Applied Computer Science, AGH

## Structure of the project:
- 'src/' - source code (main.cpp, diagonalization.hpp)
- 'data/' - solutions of the quantum harmonic oscillator equation
- 'plots/' - graphs and script for plotting the data (PNG)
- 'report/' - report of a project (PDF)

## Requirements:
- g++ with C++17 support (tested with GCC 15)
- LAPACK (header `lapacke.h`, linked with `-llapacke`);
- Ubuntu/Debian: `sudo apt install liblapack-dev`
- Python 3.8+ (tested with 3.14)
- matplotlib: `sudo apt install python3-matplotlib` or `pip install matplotlib`

## System
- OS: Ubuntu 26.04.1 LTS, Linux 7.0.0-38-generic (aarch64)
- Virtualization: QEMU virtual machine
- CPU: 4 cores, 1 thread per core, 64-bit ARM (aarch64), vendor ID reported as Apple
- RAM: 3.3 GiB (plus 3.8 GiB swap)
- Compiler: g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0
- Python: <3.14.4>
- matplotlib: <3.10.7+dfsg1>
- BLAS: OpenBLAS <0.3.32+ds-5>
- LAPACK: liblapack3 <3.12.1>

## Build and run

From the `src/` directory:

    g++ -std=c++17 -O2 -Wall -Wextra main.cpp -o main -llapacke


The program writes the norm files to `../data/`, so run it from `src/`
and make sure the `data/` directory exists.

## Plots

From the project root:

    python3 plots/plots.py --n 100 --spectrum data/spectrum.csv --vectors data/eigenvectors.csv

    python3 plots/plots.py --n 1000 --spectrum data/spectrum.csv --vectors data/eigenvectors.csv