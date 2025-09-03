# FOR_PRSA Codebase

The C++ codes in each folder implement the different steps in the algorithm
for the resummation procedure discussed in the [paper](https://doi.org/10.1098/rspa.2024.0843) for the divergent weak-
field expansion for the Heisenberg-Euler Lagrangian (HEL). Each subdirectory
contains a `README.pdf` file that gives a detailed explanation of the code.

## Directory Structure

```
FOR_PRSA/
├── Electric/
│   ├── 1_5k_mom_1_5kdig
│   ├── 1k_mom_1kdig
│   ├── 2h_mom_3kdig
│   ├── 2k_mom_2kdig
│   ├── 5h_mom_550dig
│   ├── construct_P
│   ├── Delta
│   ├── moments
│   └── Pade_1k
├── Magnetic/
│   ├── 1_5k_mom_1kdig
│   ├── 1h_mom_130dig
│   ├── 1k_mom_1kdig
│   ├── 20_mom_1hdig
│   ├── 2h_mom_250_dig
│   ├── 2k_mom_2kdig
│   ├── 50_mom_1hdig
│   ├── 5h_mom_5hdig
│   ├── 5_mom_1hdig
│   ├── construct_p
│   ├── Delta
│   ├── moments
│   └── Pade_50
```

The subdirectories `Magnetic` and `Electric` contain the source codes relevant
to the resummation of the HEL in purely magnetic and purely electric
backgrounds, respectively.

- **construct_p**: Code for constructing the matrix in equation (4.25) of the paper.  
- **Delta**: Code that implements the resummation in reference [5].  
- **moments**: Generates perturbation coefficients for the weak-field expansion in equation (3.3).  
- **Pade_50 / Pade_1k**: Implements Padé approximants of the divergent weak-field expansions for the HEL using 50 and 1000 perturbation coefficients, respectively.

The remaining directories correspond to implementations of the resummation
algorithm for various numbers of moments and working precisions.

### Example: Magnetic/1h_mom_130dig

```
Magnetic/1h_mom_130dig/
├── Constants/
├── first/
├── second/
├── third/
├── fourth/
├── function/
└── results/
```

- **Constants**: Solves the system of linear equations in equation (4.24) using 100 perturbation coefficients at 130-digit working precision.  
- **first, second, third, fourth**: Compute the respective terms in equation (4.26).  
- **function**: Computes the second term β∆(β) in equation (4.20).  

## Dependencies

This code uses the **C++ Boost Library**, specifically:

- [Boost.Multiprecision](https://www.boost.org/doc/libs/) (requires **GNU GMP**, **GNUMPFR**, **GNUMPC**)  
- [Boost.MPI](https://www.boost.org/doc/libs/) (requires an MPI implementation; tested with **OpenMPI v5.0.5**)  
- [Eigen 3](https://eigen.tuxfamily.org/) for LU factorization with `mpreal` datatypes from the [MPFR C++ library](http://www.holoborodko.com/pavel/mpfr/).  


## How to Build and Run

### Prerequisites

Make sure the following dependencies are installed on your system:

- A modern **C++ compiler** (supporting C++17 or newer)
- [Boost](https://www.boost.org/) (with Multiprecision and MPI modules)
- [GMP](https://gmplib.org/), [MPFR](https://www.mpfr.org/), and [MPC](http://www.multiprecision.org/mpc/) libraries
- [OpenMPI](https://www.open-mpi.org/) (v5.0.5 recommended)
- [Eigen 3](https://eigen.tuxfamily.org/)

On Ubuntu/Debian, you can install some prerequisites via:

```bash
sudo apt-get update
sudo apt-get install build-essential cmake libboost-all-dev libgmp-dev libmpfr-dev libmpc-dev libeigen3-dev openmpi-bin libopenmpi-dev
```

### Building

Each subdirectory contains C++ source files and a `Makefile` (or CMake configuration if provided).  
For example, to build inside `Magnetic/1h_mom_130dig`:

```bash
cd Magnetic/1h_mom_130dig
make
```

Or, if using CMake (if `CMakeLists.txt` is present):

```bash
mkdir build && cd build
cmake ..
make -j
```

### Running

Executables are usually placed in the corresponding `results/` or build folder.  
Run them with MPI as needed, for example:

```bash
mpirun -np 4 ./my_executable input_file
```

Replace `my_executable` and `input_file` with the appropriate names for the code in each subdirectory.

### Notes

- Precision settings (e.g., 130-digit working precision) are hardcoded in some files; you may adjust these in the source if required.  
- Output files are typically stored in the `results/` subfolder for each module.  
