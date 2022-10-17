//=======================================================================
//
// Example program to test out tree decomposition algorithm
//
// Reads in a path of length 10, and outputs a tree and nice tree-decomposition to separate .dot files.
//
//
//=======================================================================
#include "tree_decomposition.hpp"
#include "nice_tree_decomposition.hpp"
#include <iostream>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Graph_traits = graph_traits<Graph>;
using Graph_vertex = Graph_traits::vertex_descriptor;
using Graph_vertex_size_type = Graph_traits::vertices_size_type;
using Decomposition_bags_map = std::map<Graph_vertex, std::set<Graph_vertex>>;
using Decomposition_bags = associative_property_map<Decomposition_bags_map>;

int main() {
    Graph path = Graph(10);

    for(int i=0; i < 8; i++) {
        add_edge(i, i+1, path);
    }

    // std::ofstream path_graph;
    // path_graph.open("path.dot");
    // write_graphviz(path_graph, path);
    // path_graph.close();

    Decomposition_bags_map bm;
    Decomposition_bags b(bm);

    Graph dec_path;

    bool x = tree_decomposition(path, dec_path, b, (-0.75));

    std::cout << "this has been updated\n";

    for(int i=0; i < 10; i++) {
        std::cout << "{";
        std::set<Graph_vertex> current_bag = get(b,i);
        for (auto const &item: current_bag) {
            std::cout << item << " ";
        }
        std::cout << "}\n";
    }

    std::ofstream path_dec;
    path_dec.open("path_dec.dot");
    boost::write_graphviz(path_dec, dec_path);
    path_dec.close();


}

