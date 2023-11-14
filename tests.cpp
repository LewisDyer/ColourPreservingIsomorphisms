#define BOOST_TEST_MODULE CPITests
#include <boost/test/included/unit_test.hpp>
#include "graph_generators.h"
#include "graph_utils.h"
#include "count_colourful_isomorphisms.h"

using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;

BOOST_AUTO_TEST_CASE(test_addition) {
    int result = 2 + 3;
    BOOST_CHECK_EQUAL(result, 5);
}

BOOST_AUTO_TEST_CASE(path_1_has_0_edges) {
    Graph P =  path(1);

    BOOST_CHECK_EQUAL(boost::num_edges(P), 0);
}

BOOST_AUTO_TEST_CASE(path_3_has_2_edges) {
    Graph P =  path(3);

    BOOST_CHECK_EQUAL(boost::num_edges(P), 2);
}

BOOST_AUTO_TEST_CASE(path_100_has_99_edges) {
    Graph P =  path(100);

    BOOST_CHECK_EQUAL(boost::num_edges(P), 99);
}

BOOST_AUTO_TEST_CASE(star_1_has_0_edges) {
    Graph P =  star(1);

    BOOST_CHECK_EQUAL(boost::num_edges(P), 0);
}

BOOST_AUTO_TEST_CASE(star_3_has_2_edges) {
    Graph P =  star(3);

    BOOST_CHECK_EQUAL(boost::num_edges(P), 2);
}

BOOST_AUTO_TEST_CASE(star_100_has_99_edges) {
    Graph P =  star(100);

    BOOST_CHECK_EQUAL(boost::num_edges(P), 99);
}

BOOST_AUTO_TEST_CASE(clique_1_has_0_edges) {
    Graph P = clique(1);

    BOOST_CHECK_EQUAL(boost::num_edges(P),0);
}

BOOST_AUTO_TEST_CASE(clique_3_has_3_edges) {
    Graph P = clique(3);

    BOOST_CHECK_EQUAL(boost::num_edges(P),3);
}

BOOST_AUTO_TEST_CASE(clique_100_has_4950_edges) {
    Graph P = clique(100);

    BOOST_CHECK_EQUAL(boost::num_edges(P),4950);
}

BOOST_AUTO_TEST_CASE(single_vertex_found_in_path_1000) {
    DiGraph T = makeRooted(clique(1), 0);

    Colours col_T;

    col_T[0] = 0;

    ColourMap colour_T(col_T);

    Graph G = path(1000);

    Colours col_G;

    for(int i=0; i < boost::num_vertices(G); i++) {
        col_G[i] = 0;
    }

    ColourMap colour_G(col_G);

    BOOST_CHECK_EQUAL(new_tree_count(T, 1, colour_T, G, colour_G), 1000);
}