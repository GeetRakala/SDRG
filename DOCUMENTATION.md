# SDRG Technical Documentation

This document provides a comprehensive deep-dive into the Numerical Strong Disorder Renormalisation Group (SDRG) codebase. It details the internal data structures, algorithmic flows, and mathematical rules used to simulate random quantum systems.

## 1. Codebase Architecture

The project is structured to separate concern between graph representation, algorithmic logic, and input/output operations.

| File | Purpose |
| :--- | :--- |
| **`main.cpp`** | Entry point. Orchestrates the simulation loop: `Init -> Search -> Decimate -> Repeat`. |
| **`graph.hpp`** | Defines `BGLGraph` and `Node` structs. Wraps Boost Graph Library (BGL). |
| **`graph_init.cpp`** | Lattice generation (Chain/Square) and disorder assignment (Box/Fixed). |
| **`smart_search.cpp`** | Implements the **Activated Dijkstra** algorithm to find the strongest local scale efficiently. |
| **`smart_decimate.cpp`** | Executes the renormalization rules (Node/Edge decimation) and updates the graph. |
| **`smart_utilities.cpp`** | Helper functions for topological updates (edge reassignment, negative coupling resolution). |
| **`graph_sample.cpp`** | Probabilistic sampling for cluster statistics and entanglement entropy (Multithreaded). |
| **`dumb_search.cpp`** | Naive global search implementation for verification ($O(N)$). |
| **`dumb_decimate.cpp`** | Naive decimation logic for verification. |
| **`graph_io.cpp`** | Handles I/O: Parsing `config.txt`, writing JSON snapshots, and logging CSV results. |

---

## 2. Core Data Structures (`graph.hpp`)

The simulation relies on a custom wrapper around the Boost Graph Library (BGL).

### The `Node` Struct
Each site in the lattice is represented by a `Node` struct.
```cpp
struct Node {
  double range;       // Log-energy variable (Magnetic Moment \ln h_i or Bond Strength \ln J_{ij})
  NodeStatus status;  // Active, Inactive (decimated), or Dead
  int index;          // Unique identifier
  int clusterIndex;   // ID of the cluster this node belongs to (for tracking merged sites)
  int subsystem;      // 0 or 1 (for bipartite entanglement calculation)
};
```

### The `BGLGraph` Class
This class manages the adjacency list and ensures auxiliary maps are synchronized with the graph state.
-   **Underlying Type**: `boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Node, EdgeWeightProperty>`
-   **Optimization Maps**:
    -   `nodes_by_status`: Fast retrieval of all `Active` nodes for the search step.
    -   `nodes_by_clusterIndex`: Fast retrieval of all nodes in a cluster (used during edge decimation merging).

---

## 3. System Initialization

Before the simulation begins, the lattice creates a topology and assigns random disorder.

### Pseudocode: Initialization
```python
function InitializeLattice(L, type, disorder_type, delta):
    # 1. Create Topology
    if type == "square":
        Create LxL grid with nearest-neighbor edges
    else if type == "chain":
        Create L nodes in a line
    
    # 2. Assign Disorder
    for each node i:
        # Assign magnetic moment (h)
        h_max = CalculateHMax(delta, type)
        h = Random(0, h_max) if disorder_type == "box" else h_max
        node.range = -ln(h) # Convert to log-energy
        
    for each edge (i, j):
        # Assign bond strength (J)
        J = Random(0, 1)
        edge.weight = -ln(J) # Convert to log-energy
```

---

## 4. The Smart SDRG Algorithm

The core of the simulation is the "Smart" algorithm. Unlike naive approaches that search the global maximum energy ($O(N)$), this method searches for a **local** minimum in the renormalized energy landscape.

### High-Level Flow
```python
while NumberOfActiveNodes > 0:
    1. Search: Find the "Local Minimum" (Target for decimation).
       target = findLocalMinimum(graph)
       
    2. Check Stop Condition:
       If target.energy > max_energy_cutoff: Break
       
    3. Decimate: Apply renormalization rules.
       if target is Node:
           DecimateNode(target)
       else if target is Edge:
           DecimateEdge(target)
           
    4. Repair: Local Graph Repair (`smart_utilities.cpp`).
       Verify and fix local topology (distances, connections).
```

---

## 5. Step 1: The Search (Activated Dijkstra)

The critical optimization in this codebase is the **Activated Dijkstra** search. It finds the item with the smallest 'range' (largest energy) by exploring from a random active site.

### Pseudocode: Activated Dijkstra
```python
function ActivatedDijkstra(startNode, endNode, findToEnd):
    PriorityQueue pq
    pq.push(0, startNode)
    distances[startNode] = 0
    
    while pq is not empty:
        (dist, u) = pq.pop()
        
        if u == endNode: return (u, dist)
        
        for neighbor v of u:
            weight = edge_weight(u, v)
            if dist + weight < distances[v]:
                distances[v] = dist + weight
                pq.push(distances[v], v)
                
    return (closest_active_node, distance)

function findLocalMinimum():
    start = RandomActiveNode()
    
    # Use Dijkstra to find the nearest neighbor in the effective graph
    (next, dist1) = ActivatedDijkstra(start) 
    
    local_min = min(start.range, dist1)
    
    # If the start node is the minimum, it's a local minimum!
    if local_min == start.range: return (start, NodeDecimation)
    
    # Otherwise, traverse to the neighbor and check *its* neighborhood
    while True:
        (next_next, dist2) = ActivatedDijkstra(next)
        
        local_min = min(next.range, dist1, next_next.range, dist2)
        
        if local_min == next.range: return (next, NodeDecimation)
        if local_min == dist1: return (Edge between start and next, EdgeDecimation)
        
        # Determine direction of steepest descent
        start = next
        next = next_next
        dist1 = dist2
```

### Potential Optimizations
1.  **Temporal Locality**: The next decimation event is statistically likely to occur near the previous one (due to range reductions). Instead of restarting `findLocalMinimum` from a random node, biasing the start node selection towards the last decimated neighborhood could speed up convergence.
2.  **Parallel Independent Searches**: If the graph has multiple disjoint low-energy excitations (local minima separated by large distances), these could potentially be found and decimated in parallel, provided their "light cones" of influence do not overlap.

---

## 6. Step 2: Decimation & Renormalization

Once a target is identified, the physical degrees of freedom are integrated out (decimated), and the effective couplings are renormalized.

### Node Decimation
Occurs when a magnetic moment $\Omega = h_i$ is the largest energy scale. The site is frozen, and effective bonds are created across it.

**Renormalization Rule**:
**Renormalization Rule**:
$$
J_{new} \approx \frac{J_{ik} J_{il}}{h_i} \implies \tilde{d}_{kl} = d_{ik} + d_{il} - r_i
$$

**Pseudocode**:
```python
function DecimateNode(node_i):
    Mark node_i as Inactive
    
    For each pair of neighbors (k, l) of i:
        d_ik = distance(i, k)
        d_il = distance(i, l)
        
        # Calculate new effective distance
        new_dist = d_ik + d_il - node_i.range
        
        # Add effective edge between k and l
        AddEdge(k, l, new_dist)
        
    Remove all original edges connected to node_i
    ProcessNegativeEdges(k) # Check for stability
    UpdateAllEdgeDistancesFromNode(node_i) # Trigger local repair
```

### Edge Decimation
Occurs when a bond coupling $\Omega = J_{ij}$ is the largest energy scale. The two sites $i$ and $j$ effectively lock into a single cluster.

**Renormalization Rule**:
**Renormalization Rule**:
$$
h_{new} \approx \frac{h_i h_j}{J_{ij}} \implies \tilde{r}_{new} = r_i + r_j - d_{ij}
$$

**Pseudocode**:
```python
function DecimateEdge(node_i, node_j, dist_ij):
    # Merge clusters
    Cluster_i = node_i.clusterIndex
    Cluster_j = node_j.clusterIndex
    MergeClusters(Cluster_i, Cluster_j)
    
    # Calculate new effective moment
    new_range = node_i.range + node_j.range - dist_ij
    
    # Update node_i to represent the new cluster moment
    node_i.range = new_range
    Mark node_j as Dead
    
    # Reassign j's neighbors to i
    For each neighbor k of j:
        d_jk = distance(j, k)
        d_ik_effective = EffectiveDistance(i, k) # via Dijkstra check
        
        AddEdge(i, k, min(d_jk, d_ik_effective))
        
    Remove self-loops and duplicate edges on node_i
```

### Instability Handling (`processNegativeEdge`)
In 2D/3D, the subtraction $d_{ik} + d_{il} - r_i$ can sometimes yield a negative value, implying an unphysical "infinite" coupling. This code handles it by recursively merging nodes.

**Pseudocode**:
```python
function ProcessNegativeEdges(node):
    For each neighbor k of node:
        if distance(node, k) < 0:
            # Treats infinite coupling as an immediate edge decimation
            Make k Dead
            Merge Cluster(k) into Cluster(node)
            Recursively ProcessNegativeEdges(node)
```

---

## 7. Step 3: Local Graph Repair (`smart_utilities.cpp`)

A significant amount of development effort went into ensuring the "Smart" algorithm maintains a consistent graph state after local updates. Since we do not rebuild the graph from scratch, we must manually repair the local topology.

### 1. `updateAllEdgeDistancesFromNode`
**Problem**: When a node's moment $h_i$ changes (or is decimated), the effective distances to *all* its neighbors change.
**Solution**: This function iterates over all incident edges and updates their weights:
**Solution**: This function iterates over all incident edges and updates their weights:
$$
d_{new} = d_{old} - \frac{range_{node}}{2}
$$
It also recursively checks for negative edge creation (instabilities) and triggers `processNegativeEdge` if found.

**Pseudocode**:
```python
function UpdateAllEdgeDistancesFromNode(node):
    half_range = node.range / 2.0
    
    for each neighbor k of node:
        old_dist = distance(node, k)
        new_dist = old_dist - half_range
        
        UpdateEdge(node, k, new_dist)
        
        if new_dist < 0:
            ProcessNegativeEdge(edge(node, k))

    # Recursively process any new negative edges created
    while newly_negative_edges_exist:
        ProcessNegativeEdges(...)
```

### 2. `reassignInactiveEdges`
**Problem**: When a chain of nodes is decimated into a cluster, the original edges connecting the internal (now dead) nodes to the outside world must be rerouted to the new cluster "head".
**Solution**: The function calculates the shortest path through the decimated cluster (using Activated Dijkstra) and creates new edges from the cluster head to the external neighbors, ensuring the effective coupling is physically correct.

**Pseudocode**:
```python
function ReassignInactiveEdges(new_source, path):
    first_node = path[0]
    last_node = path.back()
    
    for k in path[1...end-1]: # Nodes inside the decimated cluster
        # Find shortest path from k to the new cluster boundary
        d1 = ActivatedDijkstra(k, first_node)
        d2 = ActivatedDijkstra(k, last_node)
        dist_to_boundary = min(d1, d2)
        
        # Move edges from k to new_source
        for each neighbor m of k:
            old_weight = distance(k, m)
            new_weight = old_weight + dist_to_boundary
            AddEdge(new_source, m, new_weight)
            RemoveEdge(k, m)
```

### 3. Topological Cleanup
-   **`removeSelfLoops`**: Merging nodes often creates edges from a node to itself ($i \to i$). These are unphysical in this model and are removed.
-   **`removeDuplicateEdges`**: Merging can create multiple parallel edges between two nodes. This function keeps only the strongest bond (minimum distance) and removes the rest.

**Pseudocode**:
```python
function RemoveSelfLoops(node):
    for each neighbor k of node:
        if k == node:
            RemoveEdge(node, k)

function RemoveDuplicateEdges(node):
    best_edge_to_neighbor = Map()
    
    # 1. Identify best edges
    for each edge e from node to k:
        if k not in best_edge_to_neighbor:
             best_edge_to_neighbor[k] = e
        else:
             if weight(e) < weight(best_edge_to_neighbor[k]):
                  best_edge_to_neighbor[k] = e
                  
    # 2. Remove all others
    for each edge e from node:
         k = target(e)
         if e != best_edge_to_neighbor[k]:
              RemoveEdge(e)
```

### 4. Potential Optimizations

Currently, the repair logic prioritizes correctness over absolute peak performance. Several algorithmic improvements could be implemented to speed up this phase:

1.  **Batch Dijkstra for `reassignInactiveEdges`**:
    -   **Current**: The code calls `ActivatedDijkstra` separately for every node inside the decimated cluster to find the distance to the boundary ($O(k \cdot D)$ where $k$ is cluster size).
    -   **Proposed**: Run Dijkstra *once* starting simultaneously from all boundary nodes neighbors inwards. This effectively computes the distance field in a single pass ($O(D)$).

2.  **Edge Deduplication via Data Structures**:
    -   **Current**: `removeDuplicateEdges` manually filters edges using a hash map after every merge operation.
    -   **Proposed**: Switch the underlying BGL graph type to use `boost::setS` for the edge container. This would automatically enforce unique edges at the data structure level, eliminating the need for manual cleanup passes (at the cost of slightly slower edge insertion).

3.  **Lazy Updates**:
    -   **Current**: `updateAllEdgeDistancesFromNode` eagerly updates every single neighbor's edge weight immediately after a renormalization step.
    -   **Proposed**: Store a "range shift" accumulator on nodes. Only apply the subtraction $d_{new} = d_{old} - \Delta range$ lazily when the edge is next accessed during a search. This avoids processing edges that might be decimated later before they are ever traversed again.

---

## 8. The Naive SDRG Algorithm ("Dumb" Mode)

For verification and performance benchmarking, the codebase includes a "Dumb" mode (enabled via `dumb = true` in `config.txt`). This implements the standard, brute-force SDRG approach.

### Search Logic (`dumb_search.cpp`)
Unlike the Smart mode which searches locally, the Naive mode scans the **entire** graph at every step to find the global minimum energy scale.
**Complexity**: $O(N)$ per step, leading to $O(N^2)$ total complexity.

### Decimation Logic (`dumb_decimate.cpp`)
The rules are identical, but no complex local repair is needed since the neighbor distances will be re-evaluated globally in the next step anyway.

**Pseudocode**:
```python
function DumbDecimate(target):
    if target is Node:
        # Simple status update, no complex local repair needed yet
        target.status = Inactive
        neighbors = GetNeighbors(target)
        for k in neighbors:
            for l in neighbors:
                 if k != l: AddEffectiveEdge(k, l)
                 
    else if target is Edge (u, v):
        # Merge u and v
        new_range = u.range + v.range - dist(u, v)
        u.range = new_range
        v.status = Dead
        
        # Naively move all neighbors of v to u
        for neighbor k of v:
            AddEdge(u, k, distance(v, k))
            
        RemoveSelfLoops(u)
```

---

## 9. Concurrency Model

The project maximizes performance using a two-tier concurrency strategy.

### Internal Threading (Shared Memory)
Used in `graph_sample.cpp` for analyzing the final state.
-   **Cluster Statistics**: The global map of clusters is partitioned into chunks. Worker threads count sizes in parallel and merge results via `std::mutex`.
-   **Entropy Calculation**: Threads process individual clusters to calculate subsystem overlaps using a thread pool pattern.

### External Parallelism (Embarrassingly Parallel)
Used for parameter sweeps via `parallel.sh`, which wraps GNU Parallel to distribute isolated simulation instances across CPU cores.

### Potential Optimizations
1.  **Lock-Free Accumulation**: `calculateClusterStatistics` uses a shared mutex to merge results from every chunk. Using thread-local accumulation maps and a single final reduction step would eliminate contention.
2.  **Work Stealing**: Static chunking does not account for variance in cluster processing time. A work-stealing queue would better balance the load across threads.
---

## 10. Visualization

The project generates artifacts that can be visualized using Python notebooks in the `notebooks/` directory.

-   **`plotter.ipynb`**: Visualizes the lattice topology and colored clusters. Requires `json=true` in config.
-   **`cluster_dist.ipynb`**: Plots statistical distributions of cluster sizes and moments from CSV logs.

---

## 11. Advanced Topics: Graph Persistence

To assist in studying the percolation properties of the random clusters, `graph_persistence.cpp` tools are available.

-   **Sampling**: The `sampleNodes` function performs a Monte Carlo sampling of the graph nodes to estimate the number of unique clusters without traversing every single node-link structure explicitly.
-   **Metric**: `computeAverageSamples` gives insight into the "cover time" or how easily the cluster structure can be probed.

**Pseudocode**:
```python
function SampleNodes(graph, num_trials):
    total_unique_clusters = 0
    
    for i = 1 to num_trials:
        # Pick random nodes
        subset = SelectRandomNodes(graph)
        
        # Count how many distinct clusters they hit
        unique_cluster_ids = Set()
        for node in subset:
            unique_cluster_ids.add(node.clusterIndex)
            
        total_unique_clusters += unique_cluster_ids.size()
        
    return total_unique_clusters / num_trials
```

### Potential Optimizations
1.  **Union-Find Data Structure**: The current sampling logic implicitly traverses the graph or relies on properties that might be slow to query. Implementing a Union-Find (Disjoint Set Union) structure would allow $O(1)$ verification of whether two nodes belong to the same cluster, significantly speeding up the unique cluster counting in `SampleNodes`.
