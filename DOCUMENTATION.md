# SDRG Technical Documentation

This document details the internal data structures, algorithmic flows, and renormalization rules implemented in the Numerical SDRG codebase.

## 1. Codebase Architecture

The project is structured to separate concerns between graph representation, algorithmic logic, and input/output operations.

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

The core of the simulation is the "Smart" algorithm. Unlike naive approaches that search the global maximum energy ( $O(N)$ ), this method searches for a **local** minimum in the renormalized energy landscape.

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

The key optimization is the **Activated Dijkstra** search: it finds the smallest 'range' (largest energy) by exploring from a random active site rather than scanning all nodes.

### `ActivatedDijkstra`
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

Once a target is identified, the corresponding degree of freedom is integrated out and effective couplings are renormalized.

### Node Decimation
Occurs when a magnetic moment $\Omega = h_i$ is the largest energy scale. The site is frozen and effective bonds are created across it.
```math
J_{new} \approx \frac{J_{ik} J_{il}}{h_i} \implies \tilde{d}_{kl} = d_{ik} + d_{il} - r_i
```
```python
function DecimateNode(node_i):
    # Mark node_i and all nodes in the same cluster as Inactive
    cluster = GetCluster(node_i)
    for node in cluster:
        Mark node as Inactive
    
    # Update all outgoing edge distances
    # This triggers renormalization of effective couplings
    UpdateAllEdgeDistancesFromNode(node_i)
```

> [!NOTE]
> In the Smart algorithm, node decimation does NOT explicitly create new edges between pairs of neighbors. Instead, the effective couplings are implicitly maintained through the `updateAllEdgeDistancesFromNode` call, which renormalizes the outgoing distances. New effective edges between former neighbors are discovered naturally by subsequent Dijkstra searches through the now-Inactive node.

### Edge Decimation
Occurs when a bond coupling $\Omega = J_{ij}$ is the largest energy scale. Sites $i$ and $j$ lock into a single cluster.
```math
h_{new} \approx \frac{h_i h_j}{J_{ij}} \implies \tilde{r}_{new} = r_i + r_j - d_{ij}
```
```python
function DecimateEdge(path, dist_ij):
    # path = [firstNode, ...intermediate inactive nodes..., lastNode]
    firstNode = path[0]
    lastNode = path.back()
    
    # 1. Calculate new effective moment for merged cluster
    new_range = firstNode.range + lastNode.range - dist_ij
    firstNode.range = new_range
    
    # 2. Merge cluster indices
    MergeClusters(firstNode.clusterIndex, lastNode.clusterIndex)
    
    # 3. Update neighbor edge distances via Dijkstra
    UpdateNeighborEdgeDistances(lastNode)
    UpdateFirstNeighborEdgeDistances(lastNode, firstNode)
    
    # 4. Reassign all edges from lastNode to firstNode
    ReassignEdges(firstNode, lastNode)
    
    # 5. Reassign edges from intermediate inactive nodes in path
    ReassignInactiveEdges(firstNode, path)
    
    # 6. Topological cleanup
    RemoveSelfLoops(firstNode)
    RemoveDuplicateEdges(firstNode)
    
    # 7. Mark all nodes except firstNode as Dead
    for node in path[1:]:
        Mark node as Dead
```

### Instability Handling (`processNegativeEdge`)
In 2D/3D, the subtraction $d_{ik} + d_{il} - r_i$ can yield a negative value (unphysical "infinite" coupling). This is resolved by creating shortcut edges that bypass the problematic node.
```python
function ProcessNegativeEdge(edge(source, target), negative_distance):
    # Create shortcut edges from source to all neighbors of target
    for each neighbor k of target (excluding source):
        d_target_k = distance(target, k)
        new_distance = negative_distance + d_target_k
        
        if edge(source, k) exists:
            # Keep the minimum distance
            UpdateEdge(source, k, min(existing_distance, new_distance))
        else:
            AddEdge(source, k, new_distance)
        
        # If new edge is also negative, queue for processing
        if new_distance < 0:
            AddToNegativeEdgeQueue(edge(source, k))
    
    # Remove all edges connected to target
    RemoveAllEdges(target)
```

> [!IMPORTANT]
> `processNegativeEdge` does NOT mark nodes as Dead or merge clusters directly. It creates shortcut edges that effectively bypass the target node, preserving the graph's connectivity while resolving the instability.

---

## 7. Step 3: Local Graph Repair (`smart_utilities.cpp`)

The "Smart" algorithm requires careful bookkeeping to maintain a consistent graph state after local updates. The following utility functions handle this.

### 1. `updateAllEdgeDistancesFromNode`
When a node's moment $h_i$ changes (or is decimated), the effective distances to all its neighbors must be updated:
```math
d_{new} = d_{old} - \frac{r_i}{2}
```
Negative results trigger `processNegativeEdge` to handle instabilities.

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
When a chain of nodes is decimated, edges from internal (now dead) nodes must be rerouted to the cluster head. The function computes shortest paths through the decimated cluster via Activated Dijkstra and creates new edges with corrected weights.

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
-   **`removeSelfLoops`**: Removes edges from a node to itself ($i \to i$), which are unphysical.
-   **`removeDuplicateEdges`**: Keeps only the minimum-distance edge when multiple parallel edges exist between two nodes.

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

The repair logic prioritizes correctness. Potential speedups:

1.  **Batch Dijkstra for `reassignInactiveEdges`**:
    -   **Current**: The code calls `ActivatedDijkstra` separately for every node inside the decimated cluster to find the distance to the boundary ( $O(k \cdot D)$ where $k$ is cluster size).
    -   **Proposed**: Run Dijkstra *once* starting simultaneously from all boundary nodes neighbors inwards. This effectively computes the distance field in a single pass ( $O(D)$ ).

2.  **Edge Deduplication via Data Structures**:
    -   **Current**: `removeDuplicateEdges` manually filters edges using a hash map after every merge operation.
    -   **Proposed**: Switch the underlying BGL graph type to use `boost::setS` for the edge container. This would automatically enforce unique edges at the data structure level, eliminating the need for manual cleanup passes (at the cost of slightly slower edge insertion).

3.  **Lazy Updates**:
    -   **Current**: `updateAllEdgeDistancesFromNode` eagerly updates every single neighbor's edge weight immediately after a renormalization step.
    -   **Proposed**: Store a "range shift" accumulator on nodes. Only apply the subtraction $d_{new} = d_{old} - \Delta r$ lazily when the edge is next accessed during a search. This avoids processing edges that might be decimated later before they are ever traversed again.

---

## 8. The Naive SDRG Algorithm ("Dumb" Mode)

For verification and performance benchmarking, the codebase includes a "Dumb" mode (enabled via `dumb = true` in `config.txt`). This implements the standard, brute-force SDRG approach.

### Search (`dumb_search.cpp`)
Scans the entire graph at every step to find the global minimum. Complexity: $O(N)$ per step → $O(N^2)$ total.

### Decimation (`dumb_decimate.cpp`)
Same renormalization rules, but no local repair needed since the next step rescans globally.
```python
function DumbDecimate(target):
    if target is Node:
        # Mark all nodes in the same cluster as Inactive
        cluster = GetCluster(target)
        for node in cluster:
            Mark node as Inactive
        
        # Create effective edges between all pairs of neighbors
        neighbors = GetNeighbors(target)
        for k in neighbors:
            for l in neighbors:
                if k != l:
                    # Calculate effective coupling through decimated node
                    new_weight = dist(target, k) + dist(target, l) - target.range
                    AddEdge(k, l, new_weight)
        
        # Remove duplicate edges on each neighbor (keep minimum)
        for k in neighbors:
            RemoveDuplicateEdges(k)
        
        # Remove all edges from decimated node
        RemoveOutEdges(target)
                 
    else if target is Edge (path from u to v):
        # Merge u and v
        new_range = u.range + v.range - dist(u, v)
        u.range = new_range
        
        # Reassign all edges from path nodes to u
        for node in path[1:]:
            ReassignEdges(u, node)
        
        # Cleanup
        RemoveSelfLoops(u)
        RemoveDuplicateEdges(u)
        
        # Mark all nodes except first as Dead
        for node in path[1:]:
            Mark node as Dead
```

---

## 9. Concurrency Model

Two-tier concurrency strategy:

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

Python notebooks in `notebooks/`:

-   **`plotter.ipynb`**: Visualizes the lattice topology and colored clusters. Requires `json=true` in config.
-   **`cluster_dist.ipynb`**: Plots statistical distributions of cluster sizes and moments from CSV logs.

---

## 11. Advanced Topics: Graph Persistence

Tools in `graph_persistence.cpp` for studying percolation properties:

-   **`sampleNodes`**: Coupon collector sampling — counts how many random node samples are needed until every cluster has been seen.
-   **`computeAverageSamples`**: Averages the above over multiple trials (expected "cover time").
```python
function SampleNodes(graph):
    all_clusters = GetUniqueClusters(graph)
    seen_clusters = Set()
    samples = 0
    nodes = GetAllNodes(graph)
    
    # Sample until all clusters are represented
    while seen_clusters.size() < all_clusters.size():
        Shuffle(nodes)  # Randomize order
        for node in nodes:
            seen_clusters.add(node.clusterIndex)
            samples += 1
            if seen_clusters.size() == all_clusters.size():
                break
    
    return samples  # Number of samples needed to hit all clusters

function ComputeAverageSamples(graph, num_trials):
    total = 0
    for i = 1 to num_trials:
        total += SampleNodes(graph)
    return total / num_trials
```

### Potential Optimizations
1.  **Union-Find Data Structure**: The current sampling logic implicitly traverses the graph or relies on properties that might be slow to query. Implementing a Union-Find (Disjoint Set Union) structure would allow $O(1)$ verification of whether two nodes belong to the same cluster, significantly speeding up the unique cluster counting in `SampleNodes`.
