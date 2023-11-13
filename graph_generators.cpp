# include <boost/graph/adjacency_list.hpp>
# include "graph_generators.h"

using namespace boost;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
/**
 * @param n, the number of vertices in the graph.
 * @return A path on n vertices.
*/
Graph path(int n) {

    Graph path(n);

    for(int i=0; i < n-1; i++) {
        add_edge(i, i+1, path);
    }

	return path;
}


/**
 * @param n, the number of vertices in the graph.
 * @return A clique on n vertices. 
*/

Graph clique(int n) {
    Graph clique(n);;
    for(int i=0; i < n; i++) {
        for (int j=0; j < i; j++) {
            add_edge(i, j, clique);
        }    
    }

    return clique;
}

/**
 * @param n, the number of vertices in the graph.
 * @return A star on n vertices, with vertex 0 as the central vertex and vertices 1 to n-1 inclusive as the leaf vertices.
*/
Graph star(int n) {
    Graph star(n);

    for (int i=1; i < n; i++) {
        add_edge(0, i, star);
    }

    return star;
}


/**
 * @param n, the number of vertices in the graph
 * @param p, the probability of each edge existing in the graph (must be between 0 and 1 inclusive)
 * 
 * @return an Erdos-Renyi graph on n vertices with each edge existing with probability p.
*/
Graph erdos_renyi(int n, double p) {
    Graph g(n);

    boost::random::mt19937 gen;
    gen.seed(std::time(0));

    boost::random::uniform_int_distribution<> dis(0, n-1);
    boost::random::uniform_real_distribution<double> real_dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < std::min(i+1+(int)(n*p), n); ++j) {
            // (i,j) is the proposed edge to add
            if (real_dis(gen) < p) {
                boost::add_edge(i, j, g);
            }
        }
    }

    return g;
}