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
# include <boost/graph/adjacency_list.hpp>
# include <boost/tuple/tuple.hpp>
# include <boost/graph/vf2_sub_graph_iso.hpp>
# include <boost/graph/subgraph.hpp>


using namespace boost;

using Graph = adjacency_list<vecS, vecS, directedS>;
using Vertex = graph_traits<Graph>::vertex_descriptor;
using Edge = graph_traits<Graph>::edge_descriptor;
using DecompositionBags = associative_property_map<std::map<Vertex, std::set<Vertex>>>;


/**
 * @param H the pattern graph.
 * @param G the data graph.
 * @param colour_H the property map containing the vertex colours of H, is a colourful k-colouring
 * @param colour_G the property map containing the vertex colours of G, is an arbitrary k-colouring
 * @return The number of subgraphs of G isomorphic to H and preserving vertex colourings
*/
template <typename Graph, typename ColourPropertyMap>
int count_colour_preserving_isomorphisms(const Graph& H, const Graph& G, const ColourPropertyMap& colour_H, const ColourPropertyMap& colour_G) {
    // get a tree decomposition of H
    // make it nice, and get its root
    // inner recursive function (see colourful_count)
    return 1;
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

/**
 * @param v the vertex we're considering in the nice tree decomposition
 * @param tree the nice tree decomposition containing v
 * @param bags the property map of bags in the tree
 * @return The type of vertex that v is in tree (either leaf, join, forget or introduce) 
*/
template <typename Graph>
VertexType getNiceType(const Vertex& v, const Graph& tree, DecompositionBags bags) {
    // Count children
    // If 0 children, return LEAF
    // If 2 children, return JOIN
    // Then compare bag of v and the bag of child of v
    // If v's bag is bigger return FORGET, else return INTRODUCE

    switch(boost::out_degree(v, tree)) {
        case 0: return leaf; break;
        case 2: return join; break;
        case 1:
            typename graph_traits<Graph>::edge_iterator child, child_end;
            boost::tie(child, child_end) = out_edges(v, tree);
            VertexType type = (size(bags[v]) > size(bags[child])) ? forget : introduce;
            return type;
            break;
    }
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
template <typename Graph, typename ColourPropertyMap, typename BagPropertyMap>
int colourful_count(const Vertex& root, const Graph& tree, std::set<Vertex> K, Graph H, Graph G, BagPropertyMap bags,  ColourPropertyMap colour_H, ColourPropertyMap colour_G ) {
    // Get the vertices in the root bag
    // Get the subgraph of H induced by these vertices

    auto root_vertices = bags[root];

    auto Xy = H.create_subgraph(boost::make_iterator_range(root_vertices.begin(), root_vertices.end()));

    // Get the subgraph of G induced by K
    auto Gk = G.create_subgraph(boost::make_iterator_range(K.begin(), K.end()));

    std::vector<int> isomorphism(boost::num_vertices(Xy));

    //auto callback = [&isomorphism](int v1, int v2) { isomorphism[v1] = v2; };
    
    bool areIsomorphic  = boost::vf2_subgraph_iso(Xy, Gk, boost::make_iterator_property_map(colour_H.begin(), boost::get(boost::vertex_index, Xy)), boost::make_iterator_property_map(colour_G.begin(), boost::get(boost::vertex_index, Gk)));

    // Check if these subgraphs are colour-preserving isomorphic
    // If not, return 0

    if (!areIsomorphic) { return 0; }

    // Check the type of the root node in the tree ("leaf"/"join"/"forget"/"introduce")

    VertexType root_type = getNiceType(root, tree, bags);

    switch(root_type) {
        case leaf:
            return 1;
            break;
        case join:
            
            int count = 0;
            Vertex child1, child2;
            for (std::pair<typename Graph::out_edge_iterator, typename Graph::out_edge_iterator> edge_pair = boost::out_edges(root, tree); edge_pair.first != edge_pair.secondl; ++edge_pair.first) {
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
        case forget:

            typename graph_traits<Graph>::edge_iterator child, child_end;
            boost::tie(child, child_end) = boost::out_edges(root, tree);
            break;
        case introduce:

            //typename graph_traits<Graph>::edge_iterator child, child_end;
            boost::tie(child, child_end) = boost::out_edges(root, tree);
            break;
    }

    // If "leaf", return 1
    // If "join", return product of counts rooted children, use the same K
    // If "forget", get the colour of the new node added in the child. Then iterate
    // over all nodes in G with the same colour, and branch here to recurse rooting on the child with all possible extensions to K.
    // If "introduce", find the vertex that was removed in the child, find its colour, then remove the vertex in K with that colour and recurse from the child.
}