# include "tree_decomposition.hpp"
# include "nice_tree_decomposition.hpp"
# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/random.hpp>
# include <boost/graph/graphviz.hpp>
# include <random>

using namespace boost;

// using Graph = adjacency_list<vecS, vecS, directedS>;
// using Vertex = graph_traits<Graph>::vertex_descriptor;
// using Edge = graph_traits<Graph>::edge_descriptor;

template <class G>
G path(int n);

template <class G, class PM>
void save_graph(std::string filename, G & g, PM & pm);

template <class G>
void save_graph(std::string filename, G & g);