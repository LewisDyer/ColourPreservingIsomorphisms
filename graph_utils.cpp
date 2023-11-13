# include <boost/graph/adjacency_list.hpp>
# include <boost/graph/graphviz.hpp>

using namespace boost;


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