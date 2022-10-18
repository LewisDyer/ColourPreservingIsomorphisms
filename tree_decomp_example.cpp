//=======================================================================
//
// Example program to test out tree decomposition algorithm
//
// Reads in a path of length 10, and outputs a tree and nice tree-decomposition to separate .dot files.
//
//
//=======================================================================
#include "tree_decomposition.hpp"
#include "nice_tree_decomposition.hpp"
#include <iostream>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Graph_traits = graph_traits<Graph>;
using Graph_vertex = Graph_traits::vertex_descriptor;
using Graph_vertex_size_type = Graph_traits::vertices_size_type;
using Decomposition_bags_map = std::map<Graph_vertex, std::set<Graph_vertex>>;
using Decomposition_bags = associative_property_map<Decomposition_bags_map>;

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

namespace std {
	template <class T>
    // Function to print out a set
	std::ostream& operator<< (std::ostream& os, const std::set<T> & in) {
		os << "[";
		for (auto it = in.begin(); it != in.end(); it++) {
            if (it != in.begin()) os << ",";
            os << *it;
        }
        os << "]";

		return os;
	}
}

void output_bags(Graph G, Decomposition_bags bags_pm) {
    for(int i=0; i < boost::num_vertices(G); i++) {
        std::cout << i << ": " << get(bags_pm,i) << "\n";
    }
}

int main() {
    Graph path = Graph(10);

    for(int i=0; i < 9; i++) {
        add_edge(i, i+1, path);
    }


    save_graph("path.dot", path);

    Decomposition_bags_map bags_m;
    Decomposition_bags bags_pm(bags_m);

    Graph dec_path;

    bool x = tree_decomposition(path, dec_path, bags_pm, (-0.75));

    save_graph("path_dec.dot", dec_path, bags_pm);

    std::cout << "this has been updated\n";

    output_bags(dec_path, bags_pm);    

    std::cout << "======\n";

    Graph nice_path;
    Decomposition_bags_map nice_bags_m;
    Decomposition_bags nice_bags_pm(nice_bags_m);

    Graph_vertex root = boost::nice_tree_decomposition(dec_path, bags_pm, nice_path, nice_bags_pm);

    save_graph("nice_path.dot", nice_path, nice_bags_pm);

    output_bags(nice_path, nice_bags_pm);

}

