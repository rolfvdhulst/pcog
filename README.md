## PCOG: Parallelized Coloring Of Graphs
This repository contains a parallelized branch-and-price code for exactly solving the graph coloring problem, also known as vertex coloring problem.

PCOG includes safe bounding techniques and heuristics for the pricing problem as described in:

Held, S., Cook, W., Sewell, E.C. (2011). Safe Lower Bounds for Graph Coloring. In: Integer Programming and Combinatoral Optimization. IPCO 2011. Lecture Notes in Computer Science, vol 6655. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-20807-2_21

Additionally, PCOG also uses and extends some of the presolving techniques described in:

Strash, D., Thompson, L. (2022). Effective Data Reduction for the Vertex Clique Cover Problem. In: 2022 Proceedings of the Symposium on Algorithm Engineering and Experiments (ALENEX), pages 41 - 53. SIAM. https://doi.org/10.1137/1.9781611977042.4

### Dependencies
The library has the following dependencies:
- [SoPlex](https://github.com/scipopt/soplex) version 7.1.3 or greater

Users who want to build the tests need to install the following package:
- [Google Test](https://github.com/google/googletest)

Users who want to build the documentation need to install the following package:
- [Doxygen](https://www.doxygen.nl/), which is commonly available through apt or other package managers: `sudo apt install doxygen`

### Building the software
1. Create a build directory

`mkdir build`

2. Navigate to the build directory

`cd build`

3. Configure the build

`cmake .. -DCMAKE_BUILD_TYPE=Release`

Optionally, users can add `-DPCOG_BUILD_TESTS=ON` to build the tests, or `-PCOG_BUILD_DOCS=ON` to build the documentation.

4. Compile:

`make`

5. (Optional) Install the library

`make install`

### Running the software

PCOG currently only accepts files that are in the DIMACS format, see https://mat.tepper.cmu.edu/COLOR/instances.html for more information.

Usage:

`pcog [inputfile] [option [value]]`

where `[inputfile]` is the path to the instance. The options and all their possible values can be listed using:
`pcog -h`.