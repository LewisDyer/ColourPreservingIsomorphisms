cmake_minimum_required(VERSION 2.8.12)
project(ColourPreservingIsomorphisms)

add_definitions("-std=c++11")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boost-treewidth)

add_executable(driver main.cpp count_colourful_isomorphisms.h graph_generators.cpp graph_utils.cpp)
target_link_libraries(driver ${CONAN_LIBS})

add_executable(tests count_colourful_isomorphisms.cpp graph_generators.cpp graph_utils.cpp tests.cpp)
target_link_libraries(tests ${CONAN_LIBS})
add_test(NAME CPITests COMMAND tests)