cmake_minimum_required(VERSION 2.8.12)
project(ColourPreservingIsomorphisms)

add_definitions("-std=c++11")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boost-treewidth)

add_executable( col-pres-iso quick_tour.cpp)
target_link_libraries(col-pres-iso ${CONAN_LIBS})