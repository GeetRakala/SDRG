# SDRG Code Documentation

This document provides a technical overview of the Numerical Strong Disorder Renormalisation Group (SDRG) codebase. It is intended for developers and researchers who wish to understand or modify the internal logic of the simulation.

## Codebase Architecture

The source code is located in the `src/` directory.

### Core Components

-   **`main.cpp`**: The entry point of the application. It handles configuration parsing, graph initialization, the main decimation loop, and result logging.
-   **`graph.hpp`**: Defines the central data structures using the Boost Graph Library (BGL). It includes the `BGLGraph` class wrapper and auxiliary structs like `Node` and `ConfigParams`.
-   **`graph_io.cpp`**: Handles Input/Output operations, including parsing `config.txt` and writing JSON/CSV files.
-   **`graph_init.cpp`**: Responsible for initializing the graph lattice (chain/square) and setting up random disorder (box/fixed distributions).

### Algorithm Implementations

The project implements two variations of the SDRG algorithm:

-   **Smart SDRG (Kovács-Iglói Algorithm)**: An optimized approach that maintains valid renormalisation steps.
    -   `smart_search.cpp`: Implements the search for the global maximum energy (smallest local minimum variable) using an activated Dijkstra approach.
    -   `smart_decimate.cpp`: Handles the decimation logic (node or edge) and graph updates (renormalizing couplings and magnetic moments) for the smart algorithm.
    -   `smart_utilities.cpp`: Helper functions for efficient graph updates (reassigning edges, removing duplicates).

-   **Dumb SDRG (Naive Algorithm)**: A simpler, brute-force approach used for testing and baseline comparisons.
    -   `dumb_search.cpp`: Scans the entire graph to find the strongest bond/moment.
    -   `dumb_decimate.cpp`: Implements the naive decimation steps.
    -   `dumb_utilities.cpp`: Helper functions for the naive approach.

### Debugging & Analysis
-   `graph_debug.cpp`: Utilities for consistency checks (e.g., finding duplicate edges).
-   `graph_sample.cpp`: Functions for calculating physical observables like cluster statistics and entanglement entropy.
-   `graph_persistence.cpp`: Utilities relating to percolation and cluster persistence.

## Data Structures (`graph.hpp`)

The core of the simulation is the `BGLGraph` class, which wraps a Boost Adjacency List.

```cpp
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,
        Node, EdgeWeightProperty> Graph;
```

### `Node` Struct
Each vertex in the graph carries the following properties:
-   `range` (`double`): Log-energy variable (magnetic moment or bond strength).
-   `status` (`NodeStatus`): `Active`, `Inactive` (decimated but part of a cluster), or `Dead`.
-   `clusterIndex` (`int`): ID of the cluster this node belongs to.
-   `subsystem` (`int`): Partition index (0 or 1) for entanglement calculations.

### `BGLGraph` Wrapper
Provides high-level methods to manipulate the graph, ensuring auxiliary maps are kept in sync:
-   `getNodesByStatus()`: Fast access to active/inactive nodes.
-   `getNodesByClusterIndex()`: Efficient retrieval of all nodes in a specific cluster.

## Algorithm Flow

### Main Loop (`main.cpp`)

1.  **Initialization**:
    -   Parse `config.txt`.
    -   Create the lattice graph (`createLatticeGraph`).
    -   Assign subsystems (`initializeSubsystem`).

2.  **Decimation Loop**:
    -   Iterate until the target energy scale is reached or the graph is fully decimated.
    -   **Search Step**: Find the "local minimum" (the strongest bond or magnetic moment).
        -   If `dumb=true`, `findMinimumEdgeOrNodeRange` scans all elements.
        -   If `dumb=false`, `findLocalMinimum` uses the optimized search.
    -   **Decimation Step**:
        -   **Node Decimation**: If a magnetic moment is the strongest scale, the site is decimated. Its effective couplings are renormalized.
        -   **Edge Decimation**: If a bond is the strongest scale, the connected sites are merged into a cluster.
    -   **Update**: The graph topology and weights are updated locally.

3.  **Analysis**:
    -   Compute cluster size distributions.
    -   Calculate entanglement entropy and mutual information based on the final cluster configurations and subsystem partitions.
    -   Write results to `csvfiles/`.

## Difference Between Smart and Dumb Modes

-   **Dumb Mode**: Naively searches for the global maximum energy at every step and performs a direct decimation. This is $O(N^2)$ or worse.
-   **Smart Mode**: Uses a local update rule and efficient data structures (inspired by Dijkstra's algorithm) to find the next decimation target without scanning the whole system. This achieves pseudo-linear time complexity $O(N \log N)$ or better for sparse graphs.

## Concurrency Model

The project uses C++11 multithreading features in `src/graph_sample.cpp` to speed up the post-processing phase.

### Cluster Statistics (`calculateClusterStatistics`)
-   **Problem**: Counting the size of thousands of clusters can be slow.
-   **Solution**: The global map of clusters is chunked into equal parts. Each thread processes a chunk to count sizes locally.
-   **Synchronization**: A `std::mutex` is used to safely merge local counts into the global `clusterSizeCount` map.

### Entanglement Entropy (`calculateEntanglementCounts`)
-   **Problem**: Calculating subsystem overlaps (Von Neumann Entropy) requires iterating through every node of every cluster.
-   **Solution**: A thread pool pattern is simulated. Threads are spawned for each cluster processing task, up to a limit defined by `config.threads`.
-   **Synchronization**: The `entanglementCounts` map is protected by a mutex to prevent race conditions when updating configuration counts.

## Advanced Algorithm Details

### Activated Dijkstra Search (`smart_search.cpp`)
Standard SDRG implementations typically search for the global maximum energy in $O(N)$ or $O(N \log N)$ using a heap. This implementation goes further by realizing that the "local minimum" in the energy landscape (conceptually equivalent to the strongest bond/moment) can be found by exploring only a localized neighborhood of the last decimated site.
-   **Method**: Uses a modified Dijkstra algorithm to explore the "active" graph.
-   **Benefit**: This reduces the effective search space significantly, as most renormalization events only affect their immediate spatial vicinity.

### Renormalization Instability Handling (`smart_utilities.cpp`)
A known issue in higher-dimensional SDRG is the generation of "negative" effective couplings during the renormalization step `new_dist = dist_1 + dist_2 - range_node`.
-   **Edge Case**: If the new distance becomes negative, it implies a breakdown of the strong disorder assumption locally or a specific topological frustration.
-   **Resolution (`processNegativeEdge`)**: The code detects these negative effective distances and immediately "processes" them by merging the connected nodes, effectively treating them as having an infinitely strong coupling (0 distance). This recursive resolution ensures the graph remains physically consistent.

### Graph Persistence (`graph_persistence.cpp`)
To assist in studying the percolation properties of the random clusters:
-   **Sampling**: The `sampleNodes` function performs a Monte Carlo sampling of the graph nodes to estimate the number of unique clusters without traversing every single node-link structure explicitly.
-   **Metric**: `computeAverageSamples` gives insight into the "cover time" or how easily the cluster structure can be probed.


