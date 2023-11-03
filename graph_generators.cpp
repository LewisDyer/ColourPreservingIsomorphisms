# include <boost/graph/adjacency_list.hpp>

# include "graph_generators.h"



using namespace boost;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
/**
 * @param n, the number of vertices in the graph.
 * @return A path on n vertices.
*/
template <class G>
G path(int n) {

    G path = Graph(n);

    for(int i=0; i < n-1; i++) {
        add_edge(i, i+1, path);
    }

	return path;
}


/**
 * @param n, the number of vertices in the graph.
 * @return A clique on n vertices. 
*/
template <class G>
G clique(int n) {
    G clique = Graph(n);
    for(int i=0; i < n; i++) {
        for (int j=0; j < i; j++) {
            add_edge(i, j, clique);
        }    
    }

    return clique;
}

template <class G>
G star(int n) {
    G star = Graph(n);

    for (int i=1; i < n; i++) {
        add_edge(0, i, star);
    }

    return star;
}

Graph erdos_renyi(int n, double p) {
    Graph g(n);

    boost::random::mt19937 gen;
    gen.seed(std::time(0));

    boost::random::uniform_int_distribution<> dis(0, n-1);
    boost::random::uniform_real_distribution<double> real_dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < std::min(i+1+(int)(n*p), n); ++j) {
            if (real_dis(gen) < p) {
                boost::add_edge(i, j, g);
            }
        }
    }

    return g;
}