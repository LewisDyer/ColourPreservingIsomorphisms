/***
 * Contains various timing methods for new_tree_count 
*/

# include "graph_utils.h"
# include "graph_generators.h"
# include "count_colourful_isomorphisms.h"

using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;