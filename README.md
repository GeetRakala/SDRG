# Numerical Strong Disorder Renormalisation Group (SDRG)

This repository contains an efficient C++ implementation of the numerical Strong Disorder Renormalisation Group (SDRG) algorithm, primarily designed for studying random quantum models such as the Ising and Potts models. The implementation leverages the Boost Graph Library for specialized graph representations and supports parallel execution for performance scaling.

## Overview

The Strong Disorder Renormalisation Group (SDRG) method is a powerful technique for studying the low-energy physics of disordered quantum many-body systems. This code implements the algorithms developed by István A. Kovács and Ferenc Iglói, as detailed in their seminal papers.

### Key Features
- **Efficient Graph Representation**: Utilizes Boost Graph Library for optimized performance.
- **Model Support**: Capable of simulating random transverse-field Ising models and Potts models.
- **Lattice Support**: Supports both **Chain** (1D) and **Square** (2D) lattice geometries.
- **Parallel Execution**: Includes scripts for parallelizing simulations across multiple parameters using GNU Parallel.
- **Advanced Metrics**: Calculates cluster statistics and entanglement measures:
  - Von Neumann Entropy
  - Mutual Information
  - Logarithmic Negativity
- **Smart & Dumb Modes**: Includes both critical (smart) and naive (dumb) decimation strategies for comparative study.

For a deep dive into the code architecture and algorithms, please refer to [DOCUMENTATION.md](DOCUMENTATION.md).

## Prerequisites

To build and run this project, you need:
- **C++ Compiler**: C++17 compliant compiler (e.g., `g++` 9.0+, `clang++`).
- **Boost C++ Libraries**: Essential for graph data structures.
  - macOS: `brew install boost`
  - Ubuntu/Debian: `sudo apt-get install libboost-all-dev`
- **GNU Parallel** (Optional): For running parallel simulations.
  - macOS: `brew install parallel`
- **Python 3**: For running visualization notebooks and scripts.

## Installation

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/GeetRakala/SDRG.git
    cd SDRG
    ```

2.  **Compile the code**:
    The project includes a `Makefile` for easy compilation.
    ```bash
    make
    ```
    This will generate an executable named `do_sdrg`.

    To clean build files:
    ```bash
    make clean
    ```

## Configuration

The simulation parameters are controlled via the `config.txt` file. You can modify this file directly or use the provided scripts to automate parameter sweeps.

| Parameter | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `dumb` | `bool` | Enable naive SDRG mode (for testing). | `false` |
| `seed` | `int` | Seed for the random number generator. | `1` |
| `graph_type` | `string` | Lattice geometry: `chain` or `square`. | `square` |
| `dist_type` | `string` | Disorder distribution: `box` or `fixed`. | `box` |
| `lattice_size` | `int` | Linear size of the lattice ($L$). Total sites $N = L$ (chain) or $L^2$ (square). | `10` |
| `delta` | `float` | Quantum control parameter (disorder strength). | `-1.0` |
| `temp` | `float` | Temperature (currently primarily for reference/cutoff). | `0` |
| `subsystems` | `int` | Number of subsystems for entanglement partitioning. | `2` |
| `persistence_trials`| `int` | Number of trials to calculate persistence of entanglement. | `1000` |
| `threads` | `int` | Number of threads for parallel cluster counting. | `10` |
| `json` | `bool` | Enable JSON output for visualization (slow). | `false` |
| `verbose` | `bool` | Enable verbose output to stdout. | `true` |

**Debug Flags**: There are several flags in `config.txt` (e.g., `debug_main`, `debug_decimate`) usually commented out, which can be enabled for granular debugging traces.

## Usage

### 1. Single Run
To run a single simulation with the current settings in `config.txt`:
```bash
./do_sdrg
# OR
make run
```
Output metrics are saved in the `csvfiles/` directory.

### 2. Parameter Sweeps (Recommended)
Use `run.sh` to execute a simulation with specific parameters in an isolated environment. This script creates a temporary directory to avoid file conflicts.

**Syntax**:
```bash
./run.sh <seed> <lattice_size> <delta>
```

**Example**:
```bash
./run.sh 42 64 0.5
```

### 3. Massively Parallel Execution
Use `parallel.sh` to run multiple instances of `run.sh` concurrently using GNU Parallel. Edit `parallel.sh` to define your parameter ranges.

```bash
# Example from parallel.sh
# Runs seeds 1 to 10 for lattice size 64 and delta 0.00
parallel ./run.sh ::: {1..10} ::: 64 ::: 0.00
```

## Project Structure

- **`src/`**: Contains all C++ source (`.cpp`) and header (`.hpp`) files.
  - `main.cpp`: Entry point and main loop.
  - `graph.hpp`: Graph definitions using Boost.
  - `smart_*.cpp`: Implementation of the optimized "Smart" SDRG algorithm.
  - `dumb_*.cpp`: Implementation of the naive "Dumb" SDRG algorithm.
- **`csvfiles/`**: Stores simulation results (Cluster statistics, Entanglement Entropy, etc.).
- **`json/`**: Stores graph state snapshots if `json=true` is enabled.
- **`scripts`**: Helper scripts (`run.sh`, `parallel.sh`, `submit.sh`).
- **`notebooks`**:
  - `plotter.ipynb`: Visualizes the SDRG graph evolution using JSON outputs.
  - `cluster_dist.ipynb`: Analyzes cluster distributions.

## Visualization & Analysis

This repository includes Python notebooks and scripts to visualize the simulation results.

### 1. Lattice Visualization (`plotter.ipynb`)

This notebook visualizes the decimation process on the 2D lattice.

-   **Requirement**: You must enable JSON output in `config.txt` by setting `json=true`.
-   **Output**: The simulation will generate `.json` files in the `json/` directory, representing the state of the graph at each step.
-   **Usage**: Open `plotter.ipynb` (which internally uses `plotter.py`) to read these JSON files and generate PDF snapshots of the lattice structure.
-   **Colors**: Clusters are color-coded to visualize how nodes are merged.

### 2. Cluster Statistics (`cluster_dist.ipynb`)

This notebook analyzes the statistical properties of the decimated graph.

-   **Input**: Reads CSV files generated in the `csvfiles/` directory (e.g., `*_statistics.csv`).
-   **Metrics**: Plots the distribution of cluster sizes and other relevant physical quantities to study the critical behavior.

## Concurrency & Performance

This project is optimized for performance using two levels of concurrency:

1.  **Internal Multithreading (Shared Memory)**:
    -   Used for calculating physical observables (Cluster Statistics, Entanglement Entropy) after the decimation is complete.
    -   Implementation: Uses `std::thread` and `std::mutex` in `src/graph_sample.cpp`.
    -   Configuration: Controlled by the `threads` parameter in `config.txt`.

2.  **External Parallelism (Embarrassingly Parallel)**:
    -   Used for running multiple independent simulations (e.g., parameter sweeps or measuring disorder averages).
    -   Implementation: Uses GNU Parallel via `parallel.sh`.
    -   Usage: `./parallel.sh` automatically distributes jobs across available CPU cores.

## Advanced Features

-   **Activated Dijkstra Optimization**: Retrieves the next decimation target using a localized shortest-path search rather than a global heap, offering superior scaling for large lattices.
-   **Instability Resolution**: Automatically detects and resolves "negative edge" instabilities (a common pathology in 2D SDRG) to maintain physical consistency.
-   **Persistence Analysis**: Includes tools to measure cluster persistence and percolation probabilities via Monte Carlo sampling.

## References

The algorithms and methodologies implemented here are based on the following research:

1.  **Paper 1**: [Phys. Rev. B 82 054437 (2010)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.054437)
2.  **Paper 2**: [Phys. Rev. B 83 174207 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.174207)
3.  **Paper 3**: [J. Phys.: Condens. Matter 23 404204 (2011)](https://iopscience.iop.org/article/10.1088/0953-8984/23/40/404204)
4.  **Paper 4**: [EPL 97 67009 (2012)](https://iopscience.iop.org/article/10.1209/0295-5075/97/67009)
5.  **Paper 5**: [Phys. Rev. B 103 174207 (2021)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.174207)

**Reviews**:
- **Part 1**: [Physics Reports 412 1 (2005)](https://www.sciencedirect.com/science/article/abs/pii/S0370157305001092)
- **Part 2**: [Eur. Phys. J. B 91 236 (2018)](https://link.springer.com/article/10.1140/epjb/e2018-90434-8)

## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.
