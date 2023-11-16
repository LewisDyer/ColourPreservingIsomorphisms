
# include <chrono>
# include "graph_utils.h"
# include "graph_generators.h"
# include "count_colourful_isomorphisms.h"
# include "performance_analysis.cpp"

using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

int main() {

    Graph uStar = star(4);
    DiGraph H = makeRooted(uStar, 1);

    Graph G = erdos_renyi(10000, 0.10);

    Colours col_H;

    // uniform colouring for H
    for(int i=0; i < boost::num_vertices(H); i++) {
        col_H[i] = i;
    }

    ColourMap colour_H(col_H);

    Colours col_G;

    // uniform colouring for G
    for(int i=0; i < boost::num_vertices(G); i++) {
    col_G[i] = i % boost::num_vertices(H);
    }

    ColourMap colour_G(col_G);

    int k = 100; //no. of executions

    for (int i=0; i < k; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        int j = new_tree_count(H, 1, colour_H, G, colour_G);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        std::cout << duration.count() << "ms \n";
    }






    return 0;
}