# ColourPreservingIsomorphisms
Implementing Lemma 1 from https://link.springer.com/chapter/10.1007/3-540-36136-7_40 in Boost.

# Installation instructions

Set up a build directory and cd into it

`mkdir build && cd build`

Then install dependencies for the project (note I don't think you need the `build` arguments here, but this is quick and dirty for now...)
`conan install .. --build=boost --build=bzip2 --build=zlib`

Configure the CMake project:

`cmake ..`

Build the application:

`cmake --build . --config Release`

