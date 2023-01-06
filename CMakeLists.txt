cmake_minimum_required(VERSION 2.8.12)
project(ColourPreservingIsomorphisms)

add_definitions("-std=c++11")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boost-treewidth)

# add_executable(quick_tour quick_tour.cpp)
# target_link_libraries(quick_tour ${CONAN_LIBS})

#add_executable(tree_example tree_decomp_example.cpp tree_decomposition.hpp nice_tree_decomposition.hpp)
#target_link_libraries(tree_example ${CONAN_LIBS})

add_executable(col_pres_iso count_colourful_isomorphisms.cpp)
target_link_libraries(col_pres_iso ${CONAN_LIBS})