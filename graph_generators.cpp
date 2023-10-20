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