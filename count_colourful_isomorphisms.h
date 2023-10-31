
# include <boost/graph/adjacency_list.hpp>
# include <boost/tuple/tuple.hpp>
# include <boost/graph/isomorphism.hpp>
# include <boost/graph/graphviz.hpp>
# include <boost/graph/undirected_graph.hpp>
# include <boost/graph/subgraph.hpp>
# include <boost/graph/copy.hpp>

using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

template <class G>
G star(int n);

template <class G>
G path(int n);

Graph erdos_renyi(int n, double p);

DiGraph makeRooted(const Graph& G, Vertex root);

template <class G, class PM>
void save_graph(std::string filename, G & g, PM & pm);

template <class G>
void save_graph(std::string filename, G & g);

int count_tree(Vertex root, Vertex rootMapped, DiGraph tree, ColourMap colour_H, Graph G, ColourMap colour_G);

int tree_count(DiGraph tree, Vertex root, ColourMap colour_H, Graph G, ColourMap colour_G);