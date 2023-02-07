# include "tree_decomposition.hpp"
# include "nice_tree_decomposition.hpp"
# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/random_spanning_tree.hpp>
# include <boost/graph/graphviz.hpp>
# include <random>


using namespace boost;

using Graph = adjacency_list<vecS, vecS, directedS>;
using Vertex = graph_traits<Graph>::vertex_descriptor;
using Edge = graph_traits<Graph>::edge_descriptor;

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

template <class G, class PM>
// Helper function to save graph to file from https://github.com/Cynt3r/boost-treewidth
void save_graph(std::string filename, G & g, PM & pm) {
	std::ofstream file;
  	file.open(filename);
  	boost::write_graphviz(file, g, make_label_writer(pm));
  	file.close();
}

template <class G>
// Same helper function but with default labels
void save_graph(std::string filename, G & g) {
	std::ofstream file;
  	file.open(filename);
  	boost::write_graphviz(file, g);
  	file.close();
}