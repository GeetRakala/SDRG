# Numerical Strong Disorder Renormalisation Group (SDRG)

C++ implementation of the numerical Strong Disorder Renormalisation Group (SDRG) algorithm for studying random quantum Ising and Potts models. Uses Boost Graph Library for graph representation and supports parallel execution.

## Overview

The SDRG method is a real-space renormalization technique for disordered quantum many-body systems. This implementation follows the algorithms of Kovács and Iglói.

### Key Features
- **Graph Representation**: Boost Graph Library (BGL) adjacency list.
- **Model Support**: Capable of simulating random transverse-field Ising models and Potts models.
- **Lattice Support**: Supports both **Chain** (1D) and **Square** (2D) lattice geometries.
- **Parallel Execution**: Includes scripts for parallelizing simulations across multiple parameters using GNU Parallel.
- **Advanced Metrics**: Calculates cluster statistics and entanglement measures:
  - Von Neumann Entropy
  - Mutual Information
  - Logarithmic Negativity
- **Smart & Dumb Modes**: Optimized local search (smart) and naive global search (dumb) for benchmarking.

See [DOCUMENTATION.md](DOCUMENTATION.md) for algorithmic details.

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

2.  **Compile**:
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
- **`scripts/`**: Helper scripts (`run.sh`, `parallel.sh`, `submit.sh`).
- **`notebooks`**:
  - `plotter.ipynb`: Visualizes the SDRG graph evolution using JSON outputs.
  - `cluster_dist.ipynb`: Analyzes cluster distributions.

## Visualization & Analysis

Python notebooks for visualization:

### 1. Lattice Visualization (`plotter.ipynb`)

Visualizes the decimation process on the 2D lattice.

-   **Requirement**: Enable JSON output (`json=true` in `config.txt`).
-   **Output**: `.json` files in `json/` representing graph state at each step.
-   **Usage**: Run `plotter.ipynb` to generate PDF snapshots.
-   **Colors**: Clusters are color-coded.

### 2. Cluster Statistics (`cluster_dist.ipynb`)

Analyzes statistical properties of the decimated graph.

-   **Input**: CSV files from `csvfiles/` (e.g., `*_statistics.csv`).
-   **Metrics**: Cluster size distributions and critical behavior.

## Concurrency & Performance

Two levels of concurrency:

1.  **Internal Multithreading (Shared Memory)**:
    -   Used for calculating physical observables (Cluster Statistics, Entanglement Entropy) after the decimation is complete.
    -   Implementation: Uses `std::thread` and `std::mutex` in `src/graph_sample.cpp`.
    -   Configuration: Controlled by the `threads` parameter in `config.txt`.

2.  **External Parallelism (Embarrassingly Parallel)**:
    -   Used for running multiple independent simulations (e.g., parameter sweeps or measuring disorder averages).
    -   Implementation: Uses GNU Parallel via `parallel.sh`.
    -   Usage: `./parallel.sh` automatically distributes jobs across available CPU cores.

## Advanced Features

-   **Activated Dijkstra**: Finds next decimation target via localized shortest-path search instead of global scan.
-   **Instability Resolution**: Detects and resolves negative-edge instabilities (2D SDRG pathology).
-   **Persistence Analysis**: Monte Carlo sampling for cluster persistence and percolation probabilities.

## References

Based on:

1.  [Phys. Rev. B 82, 054437 (2010)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.054437)
2.  [Phys. Rev. B 83, 174207 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.174207)
3.  [J. Phys.: Condens. Matter 23, 404204 (2011)](https://iopscience.iop.org/article/10.1088/0953-8984/23/40/404204)
4.  [EPL 97, 67009 (2012)](https://iopscience.iop.org/article/10.1209/0295-5075/97/67009)
5.  [Phys. Rev. B 103, 174207 (2021)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.174207)

Reviews:
-   [Physics Reports 412, 1 (2005)](https://www.sciencedirect.com/science/article/abs/pii/S0370157305001092)
-   [Eur. Phys. J. B 91, 236 (2018)](https://link.springer.com/article/10.1140/epjb/e2018-90434-8)

## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.
