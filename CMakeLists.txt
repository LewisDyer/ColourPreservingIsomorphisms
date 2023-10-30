cmake_minimum_required(VERSION 2.8.12)
project(ColourPreservingIsomorphisms)

add_definitions("-std=c++11")

# include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
# )
# # For Windows: Prevent overriding the parent project's compiler/linker settings
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# FetchContent_MakeAvailable(googletest)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boost-treewidth)



#add_executable(quick_tour quick_tour.cpp)
#target_link_libraries(quick_tour ${CONAN_LIBS})

#add_executable(tree_example tree_decomp_example.cpp tree_decomposition.hpp nice_tree_decomposition.hpp)
#target_link_libraries(tree_example ${CONAN_LIBS})

add_executable(col_pres_iso count_colourful_isomorphisms.cpp tree_decomposition.hpp nice_tree_decomposition.hpp graph_generators.cpp)
target_link_libraries(col_pres_iso ${CONAN_LIBS})

# enable_testing()

# add_executable(
#   hello_test
#   tests/tests.cc
# )
# target_link_libraries(
#   hello_test
#   GTest::gtest_main
# )

# include(GoogleTest)
# gtest_discover_tests(hello_test)