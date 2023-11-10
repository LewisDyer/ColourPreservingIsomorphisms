# ColourPreservingIsomorphisms
Given a tree T with a colourful k-vertex colouring, and a graph G with a k-vertex colouring (not necessarily colourful), returns the number of colour-preserving subgraph isomorphsms from T into G.

# Project structure

* `graph_generators` contains several methods for creating various common graph classes, such as paths or Erdos-Renyi graphs.
* `graph_utils` contains several useful methods for working with graphs, such as outputting graphs to `.dot` files.
* `count_colourful_isomorphisms` contains the primary methods for counting colour-preserving subgraph isomorphisms.
* `main.cpp` is a driver program, showing an example of creating graphs and vertex colourings and counting colour-preserving subgraph isomorphisms.
* `tests.cpp` contains unit tests.

Note that `tree_decomposition` and `nice_tree_decomposition` are used to compute tree decompositions and nice tree decompositions of pattern graphs respectively, as implemented by Václav Král at https://github.com/Cynt3r/boost-treewidth. However these are currently unused on this branch of the program.

# Installation instructions

Set up a build directory and cd into it: `mkdir build && cd build`

Then install dependencies for the project: `conan install .. --build=boost --build=bzip2 --build=zlib`

Configure the CMake project: `cmake ..`

Build the application: `cmake --build . --config Release`

