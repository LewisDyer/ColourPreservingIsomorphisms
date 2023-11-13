#ifndef GRAPH_UTILS
#define GRAPH_UTILS

#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace boost;

using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;

DiGraph makeRooted(const Graph& G, Vertex root);

template <class G, class PM>
void save_graph(std::string filename, G & g, PM & pm) {
    std::ofstream file;
    file.open(filename);
    boost::write_graphviz(file, g, make_label_writer(pm));
    file.close();
}

template <class G>
void save_graph(std::string filename, G & g) {
    std::ofstream file;
    file.open(filename);
    boost::write_graphviz(file, g);
    file.close();
}

#endif