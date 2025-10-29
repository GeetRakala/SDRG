# Numerical Strong Disorder Renormalisation Group algorithm (SDRG)

This is an efficient C++ implementation of the numerical Strong Disorder Renormalisation Group algorithm for random models (Ising/Potts) developed by István A. Kovács and Ferenc Iglói
documented in the following papers:
Paper 1: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.054437
Paper 2: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.174207
Paper 3: https://iopscience.iop.org/article/10.1088/0953-8984/23/40/404204
Paper 4: https://iopscience.iop.org/article/10.1209/0295-5075/97/67009
Paper 5: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.174207

A review of the strong disorder renormalisation method can be found here:
Part 1: https://www.sciencedirect.com/science/article/abs/pii/S0370157305001092
Part 2: https://link.springer.com/article/10.1140/epjb/e2018-90434-8

This repository uses templates provided by Boost Graph Library for efficient graph representations. 
Some bits have been parallelised using concurrency in a thread-safe manner.
Start by looking at src/main.cpp to understand the flow of the algorithm.

Repository stucture:
1. '/src' contains all the .cpp and .hpp files.
2. Compilation is handled with the 'Makefile'.
3. '/build' holds temporary build files. Can be cleared using 'make clean'.
4. The 'config.txt' file provides neccesary parameters.
5. A JSON file with the state of the graph can be printed into a 'graph_data.json' file.
6. JSON files in /json can be visualised with  'plotter.ipynb'
7. 'run.sh' isolates the executable and config files into a new unique folder.
8. 'parallel.sh' parallelises 'run.sh' using GNU Parallel.

Enjoy!
