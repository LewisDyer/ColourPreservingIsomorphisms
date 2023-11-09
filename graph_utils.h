// using Graph = adjacency_list<vecS, vecS, directedS>;
// using Vertex = graph_traits<Graph>::vertex_descriptor;
// using Edge = graph_traits<Graph>::edge_descriptor;

template <class G, class PM>
void save_graph(std::string filename, G & g, PM & pm);

template <class G>
void save_graph(std::string filename, G & g);