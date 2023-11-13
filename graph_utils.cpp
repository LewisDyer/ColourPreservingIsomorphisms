# include <string>
# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/graphviz.hpp>
# include <queue>

using namespace boost;

using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;

DiGraph makeRooted(const Graph& G, Vertex root) {
    DiGraph d;
    std::vector<bool> visited(boost::num_vertices(G), false);
    std::queue<int> q;

    visited[root] = true;
    q.push(root);

    while (!q.empty()) {
        int currentVertex = q.front();
        q.pop();

        Graph::adjacency_iterator vi, vi_end;
        for (boost::tie(vi, vi_end) = boost::adjacent_vertices(currentVertex, G); vi != vi_end; ++vi) {
            if (!visited[*vi]) {
                visited[*vi] = true;
                q.push(*vi);
                boost::add_edge(currentVertex, *vi, d);
            }
        }
    }

    return d;
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
// Same helper function, but with default labels
void save_graph(std::string filename, G & g) {
	std::ofstream file;
  	file.open(filename);
  	boost::write_graphviz(file, g);
  	file.close();
}