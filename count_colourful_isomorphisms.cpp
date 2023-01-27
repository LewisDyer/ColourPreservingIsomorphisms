// GENERAL PIPELINE
// Read in the pattern graph H, and generate a nice tree decomposition, returning the decomposition and its associated bags
// Find the root arbitrarily, get the colours of vertices in H at the root bag, then get a list of all subsets of vertices in G of the same size with the same set of colours.
// For each such subset:
//      Check if there's a subgraph isomorphism with the bag at the current node. If there isn't, just return 0.
//      Then consider the type of the current node in the nice tree decomposition:
//          LEAF: Just return 1!
//          JOIN: Take the product of counts of the two children, with the same K!
//          FORGET: Get the colour of the new node added in H, then sum together counts of that child when adding a new vertex to K with the same colour.
//          INTRODUCE: Get the colour of the node removed in H, then remove the vertex in K with the same colour and recurse.
//      Use this recursion to tot up the total number of extensions at the root.
//      If you ever get a value of at least 1, terminate iterating through Ks and return true. Else if there's no more Ks to check, return False

// BOOST-SPECIFIC ISSUES
// How do we store graph colours? I'm guessing some sort of property map
// How do we check graph isomorphisms? Keeping in mind I need to store the map that defines such an isomorphism as well for forget and introduce nodes.
# include <set>
# include <list>
# include <cmath>
# include <boost/graph/adjacency_list.hpp>
# include <boost/tuple/tuple.hpp>
# include <boost/graph/isomorphism.hpp>
# include <boost/graph/graphviz.hpp>
# include <iostream>

# include "tree_decomposition.hpp"
# include "nice_tree_decomposition.hpp"


using namespace boost;

using Graph = adjacency_list<listS, vecS, directedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex = graph_traits<Graph>::vertex_descriptor;
using Edge = graph_traits<Graph>::edge_descriptor;
using BagsMap = std::map<Vertex, std::set<Vertex>>;
using DecompositionBags = associative_property_map<BagsMap>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

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

void output_bags(Graph G, DecompositionBags bags_pm) {
    for(int i=0; i < boost::num_vertices(G); i++) {
        std::cout << i << ": " << bags_pm[i] << "\n";
    }
}

/**
 * @param n, the number of vertices in the graph.
 * @return A path on n vertices.
*/
Graph path(int n) {

    Graph path = Graph(n);

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
  	boost::write_graphviz(file, g, boost::make_label_writer(pm));
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

/**
 * @param tree a Graph that's a tree
 * @return the vertex descriptor of the root of the tree.
 * 
 * Note default behaviour for graphs with no sources is unimplemented.
 * 
*/
template <typename Graph>
Vertex get_root(const Graph& tree) {
    typedef typename graph_traits<Graph>::vertex_iterator iter_v;
    for (std::pair<iter_v, iter_v> p = vertices(tree); p.first != p.second; ++p.first) {
        if (in_degree(p.first, tree) == 0) {
            return p.first;
        }
    }
}



enum VertexType {leaf, join, forget, introduce};
const char* node_type[4] = {"leaf", "join", "forget", "introduce"};

/**
 * @param v the vertex we're considering in the nice tree decomposition
 * @param tree the nice tree decomposition containing v
 * @param bags the property map of bags in the tree
 * @return The type of vertex that v is in tree (either leaf, join, forget or introduce) 
*/
VertexType getNiceType(const Vertex& v, const Graph& tree, DecompositionBags bags) {
    // Count children
    // If 0 children, return LEAF
    // If 2 children, return JOIN
    // Then compare bag of v and the bag of child of v
    // If v's bag is bigger return FORGET, else return INTRODUCE

    switch(boost::out_degree(v, tree)) {
        case 0: { return leaf; break; }
        case 2: { return join; break; }
        case 1: {
            typename graph_traits<Graph>::out_edge_iterator child, child_end;
            boost::tie(child, child_end) = out_edges(v, tree);
            VertexType type = (size(bags[v]) > size(bags[boost::source(*child, tree)])) ? introduce : forget;
            return type;
            break;
        }
        default: { return leaf; break; }
    }
}

/**
 * @param v1 a vertex in G
 * @param v2 a vertex in H
 * @param colour1 a property map of colours in G
 * @param colour2 a property map of colours in H
 * @return true if v1 and v2 have the same colours in their respective graphs, and false otherwise.
*/
bool colour_match(Vertex v1, Vertex v2, ColourMap colour1, ColourMap colour2) {
    return colour1[v1] == colour2[v2];
}

Graph induced_subgraph(Graph G, std::set<Vertex> V) {
    Graph S(V.size());

    for (Vertex v: V) {
        Vertex v_sub = vertex(v, S);
        auto adj_vertices = boost::adjacent_vertices(v, G);
        for (auto v_target: boost::make_iterator_range(adj_vertices.first, adj_vertices.second)) {

            if (std::find(V.begin(), V.end(), v_target) != V.end()) {
                add_edge(v, v_target, S);
            }
            
        }
    }

    return S;
}

/**
 * @param root A vertex of the root of the subtree of the nice tree decomposition we're considering. Its bag X_y has m vertices, with H_y being the subgraphs of H induced by all vertices in this subtree.
 * @param tree the nice tree decomposition we're working with.
 * @param K a subset of the vertex set of G with exactly m vertices.
 * @param H a pattern graph with k vertices
 * @param G a data graph
 * @param bags a property map containing the bags of each vertex in tree
 * @param colour_H a colourful k-colouring of the vertices of H
 * @param colour_G an arbitrary k-colouring of vertices of G.
 * @return The number of colour-preserving isomorphisms between subgraphs of G and H_y that extend colour-preserving isomorphisms between K and X_y. If no such isomorphisms exist between K and X_y, return 0.
*/
int colourful_count(Vertex root, Graph tree, std::set<Vertex> K, Graph H, Graph G, DecompositionBags bags,  ColourMap colour_H, ColourMap colour_G ) {
    // Get the vertices in the root bag
    // Get the subgraph of H induced by these vertices

    std::set<Vertex> root_vertices = bags[root];

    std::cout << "root is " << root_vertices << "\n";

    std::cout << "K is " << K << "\n";


    Graph Xy = induced_subgraph(H, root_vertices);

    Graph Gk = induced_subgraph(G, K);

    save_graph("Test_Xy.dot", Xy);

    save_graph("Test_Gk.dot", Gk);

    std::vector<Vertex> iso(num_vertices(Xy));

    std::cout << "before isomorphism check\n";

    bool areIsomorphic;

    if ((num_vertices(Xy) == 0) && (num_vertices(Gk) == 0)) {
        areIsomorphic = true;
    } else {
        areIsomorphic = boost::isomorphism(Xy, Gk, isomorphism_map(make_iterator_property_map(iso.begin(), get(vertex_index, Xy), iso[0])));
    }
  
    //bool areIsomorphic = isomorphism(Xy, Gk, isomorphism_map(make_iterator_property_map(isomorphism.begin(), get(vertex_index, Xy), isomorphism[0])));
    // bool areIsomorphic = false;

    std::cout << "after iso check\n" << areIsomorphic << "\n";

    if (!areIsomorphic) { return 0; }

    // // Check the type of the root node in the tree ("leaf"/"join"/"forget"/"introduce")

     VertexType root_type = getNiceType(root, tree, bags);

     std::cout << "root type is " << node_type[root_type] << "\n" ;

     switch(root_type) {
        case leaf: {
            return 1;
            break;
        }
        case join: {
            
            int count = 0;
            Vertex child1, child2;
            for (std::pair<typename Graph::out_edge_iterator, typename Graph::out_edge_iterator> edge_pair = boost::out_edges(root, tree); edge_pair.first != edge_pair.second; ++edge_pair.first) {
                Edge e = *edge_pair.first;
                if (count == 0) {
                    child1 = boost::target(e, tree);
                } else {
                    child2 = boost::target(e,tree);
                }

                count++;
            }
            return (colourful_count(child1, tree, K, H, G, bags, colour_H, colour_G) * colourful_count(child2, tree, K, H, G, bags, colour_H, colour_G));
            break;
        }
        case forget: {

            typename graph_traits<Graph>::out_edge_iterator child, child_end;
            boost::tie(child, child_end) = boost::out_edges(root, tree);

            Vertex child1 = boost::target(*child, tree);

            //Get the colour of the new vertex in this node
            std::set<Vertex> diff;
            std::set_difference(get(bags, child1).begin(), get(bags, child1).end(), get(bags, root).begin(), get(bags, root).end(), std::inserter(diff, diff.begin()));
            Vertex new_v = *(diff.begin());
            int colour = colour_H[new_v];

            // Iterate over vertices in G and check their colour - branch on all vertices with the same colour as the new vertex.

            int total = 0;
            typedef typename graph_traits<Graph>::vertex_iterator iter_v;
            for (std::pair<iter_v, iter_v> p = vertices(tree); p.first != p.second; ++p.first) {
                if (colour_G[*p.first] == colour) {
                    std::set<Vertex> new_K;
                    std::copy(K.begin(), K.end(), std::inserter(new_K, new_K.begin()));
                    new_K.insert(*p.first);

                    total += colourful_count(child1, tree, new_K, H, G, bags, colour_H, colour_G);
                }
            }

            return total;
            break;
        }
        case introduce: {

            typename graph_traits<Graph>::out_edge_iterator child, child_end;
            boost::tie(child, child_end) = boost::out_edges(root, tree);

            Vertex child1 = boost::target(*child, tree);

            std::set<Vertex> diff;
            std::set_difference(get(bags, root).begin(), get(bags, root).end(), get(bags, child1).begin(), get(bags, child1).end(), std::inserter(diff, diff.begin()));
            Vertex forgotten_v = *(diff.begin());
            int prev_target = find(iso.begin(), iso.end(), forgotten_v) - iso.begin();

            std::set<Vertex> new_K;
            std::copy(K.begin(), K.end(), std::inserter(new_K, new_K.begin()));
            new_K.erase(prev_target);
            return colourful_count(child1, tree, new_K, H, G, bags, colour_H, colour_G);
            break;
        }
        default: {
            return 0;
            break;
        }
    }

    // If "leaf", return 1
    // If "join", return product of counts rooted children, use the same K
    // If "forget", get the colour of the new node added in the child. Then iterate
    // over all nodes in G with the same colour, and branch here to recurse rooting on the child with all possible extensions to K.
    // If "introduce", find the vertex that was removed in the child, find its colour, then remove the vertex in K with that colour and recurse from the child.


    return 1;
}

/**
 * @param H the pattern graph.
 * @param G the data graph.
 * @param colour_H the property map containing the vertex colours of H, is a colourful k-colouring
 * @param colour_G the property map containing the vertex colours of G, is an arbitrary k-colouring
 * @return The number of subgraphs of G isomorphic to H and preserving vertex colourings
*/
int count_colour_preserving_isomorphisms(Graph H, Graph G, ColourMap colour_H, ColourMap colour_G) {
    // get a tree decomposition of H
    // make it nice, and get its root
    // inner recursive function (see colourful_count)

    Graph dec;
    BagsMap bags_m;
    DecompositionBags bags_pm(bags_m);
    
    bool x = boost::tree_decomposition(H, dec, bags_pm, -0.75);

    std::cout << "initial decomposition fine\n";

    Graph tree;
    BagsMap nice_bags_m;
    DecompositionBags bags(nice_bags_m);

    Vertex root = boost::nice_tree_decomposition(dec, bags_pm, tree, bags);

    //root always has empty bags, so K starts empty
    std::set<Vertex> K;
    return colourful_count(root, tree, K, H, G, bags, colour_H, colour_G);
    //return 1;
}


int main() {

Graph H = path(5);
Graph G = path(10);

Colours col_H;

for(int i=0; i < boost::num_vertices(H); i++) {
    col_H[i] = i;
}

ColourMap colour_H(col_H);

Colours col_G;

for(int i=0; i < boost::num_vertices(G); i++) {
    if (i <=4) {
        col_G[i] = i;
    } else{
        col_G[i] = 5 - std::abs(i-4);
    }
   
}

ColourMap colour_G(col_G);

save_graph("H_colour.dot", H, colour_H);

save_graph("G_colour.dot", G, colour_G);

int count = count_colour_preserving_isomorphisms(H, G, colour_H, colour_G);

std::cout << count;

}