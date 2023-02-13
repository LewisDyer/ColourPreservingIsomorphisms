# include "graph_generators.h"

/**
 * @param n, the number of vertices in the graph.
 * @return A path on n vertices.
*/
template <class G>
G path(int n) {

    G path = G(n);

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
    G clique = G(n);
    for(int i=0; i < n; i++) {
        for (int j=0; j < i; j++) {
            add_edge(i, j, clique);
        }    
    }

    return clique;
}