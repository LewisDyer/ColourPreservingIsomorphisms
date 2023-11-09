#ifndef GRAPH_GENS
#define GRAPH_GENS

# include <set>
# include <list>
# include <queue>
# include <cmath>
# include <boost/graph/adjacency_list.hpp>
# include <boost/tuple/tuple.hpp>
# include <boost/graph/isomorphism.hpp>
# include <boost/graph/graphviz.hpp>
# include <boost/graph/undirected_graph.hpp>
# include <boost/graph/subgraph.hpp>
# include <boost/graph/copy.hpp>
# include <boost/bind.hpp>  
# include <boost/random.hpp>
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>
# include <boost/graph/breadth_first_search.hpp>
# include <iostream>

using namespace boost;

using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;

Graph path(int n);

Graph star(int n);

Graph clique(int n);

Graph erdos_renyi(int n, double p);

#endif