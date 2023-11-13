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
# include "graph_utils.h"
# include "graph_generators.h"



using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

// int count_tree(Vertex root, Vertex rootMapped, DiGraph tree, ColourMap colour_H, Graph G, ColourMap colour_G) {

//     if (boost::out_degree(root, tree) == 0 ) { //leaf node
//         return 1;
//     }
    
//     int combos = 1;
//     // iterate through children of r

//     typename graph_traits<DiGraph>::out_edge_iterator child, child_end;
//     typename graph_traits<DiGraph>::out_edge_iterator target, target_end;

//     for(boost::tie(child, child_end) = boost::out_edges(root, tree); child != child_end; ++child) {
//         Vertex childVertex = boost::target(*child, tree);
//         int childDegree = boost::degree(childVertex, tree);
//         int childColour = colour_H[childVertex];
//         int currentTotal = 0;
//         //iterate through adjacent vertices of rootMapped, look for vertices with childColour
        
//         std::pair<Graph::adjacency_iterator, Graph::adjacency_iterator> neighbors = boost::adjacent_vertices(rootMapped, G);

//         // Iterate through adjacent vertices
//         for (Graph::adjacency_iterator neighbor = neighbors.first; neighbor != neighbors.second; ++neighbor) {
//             Graph::vertex_descriptor adjacent_vertex = *neighbor;
//             int targetColour = colour_G[adjacent_vertex];
//             if (targetColour == childColour) {
//                 int targetDegree = boost::degree(adjacent_vertex, G);
//                 if (targetDegree >= childDegree) { // degree filtering to avoid checking impossible matches
//                 currentTotal += count_tree(childVertex, adjacent_vertex, tree, colour_H, G, colour_G);
//                 } else {
//                     std::cout << "avoided due to degree filtering \n";
//                 }
//             }
//         }

//         if (currentTotal == 0) {return 0;}

//         combos *= currentTotal;

//     }
    
//     return combos;

// }

// /**
//  * @param tree the pattern graph being considered (must be a tree on k vertices)
//  * @param root The root of the tree
//  * @param colour_H the colourful vertex k-colouring of the vertices of the tree
//  * @param G the data graph
//  * @param colour_G an arbitrary k-colouring of vertices of G
//  * @return The number of colour-preserving isomorphisms between H and subgraphs of G.
// */
// int tree_count(DiGraph tree, Vertex root, ColourMap colour_H, Graph G, ColourMap colour_G) {

//     //Vertex root = get_root(tree);
//     int targetColour = colour_H[root];
//     int rootDegree = boost::degree(root, tree);

//     // search for targetColour in G to consider starting points

//      int total = 0;
//             typedef typename graph_traits<Graph>::vertex_iterator iter_v;
//             for (std::pair<iter_v, iter_v> p = vertices(G); p.first != p.second; ++p.first) {
//                 if (colour_G[*p.first] == targetColour && *p.first != boost::num_vertices(G)) {
//                     if (boost::degree(*p.first, G) >= rootDegree) {
//                     total += count_tree(root, *p.first, tree, colour_H, G, colour_G);
//                 }
//                 }
//             }


//     return total;
// }

/**
 * @param tree, a directed graph of the tree pattern being considered such that there exists a path from the root to every other vertex.
 * @param v the vertex whose height is currently being calculated
 * @param heights a lookup table containing the height of each vertex so far
*/
void calculateHeights(const DiGraph& tree, Vertex v, std::vector<int>& heights) {
    int height = 0; // Initialize the height for the current vertex, initially 0 for leaf vertices
    
    // Iterate through all adjacent vertices
    for (auto it = boost::adjacent_vertices(v, tree); it.first != it.second; ++it.first) {
        int adjacent_vertex = *it.first;
        calculateHeights(tree, adjacent_vertex, heights);
        height = std::max(height, heights[adjacent_vertex] + 1); // height of a vertex is the maximum height of its children, plus one
    }

    heights[v] = height; // Set the height of the current vertex
}

/**
 * @param tree a directed graph of the tree pattern being considered such that there exists a path from the root to every other vertex.
 * @param root the root vertex of tree
 * @return A vector containing the height of every vertex in the tree.
*/
std::vector<int> getAllHeights(const DiGraph& tree, Vertex root) {
    // wrapper function for calculateHeights recursive function
    std::vector<int> heights(boost::num_vertices(tree), 0);
    calculateHeights(tree, root, heights);
    return heights;
}
/**
 * @param tree, the tree pattern being searched for
 * @param root the root vertex of the tree pattern
 * @param colour_H a property map defining a colourful k-vertex colouring of the tree
 * @param G an undirected data graph
 * @param colour_G a property map defining a k-vertex colouring of G (not necessarily colourful)
 * @return The number of colour-preserving graph isomorphisms from tree to subgraphs of G.
*/
int new_tree_count(DiGraph tree, Vertex root, ColourMap colour_tree, Graph G, ColourMap colour_G) {
    int total = 0;
    Colours cpi_count; //stores counts of partial solutions for all vertices in G

    typedef typename graph_traits<Graph>::vertex_iterator iter_v;
    std::vector<int> heights = getAllHeights(tree, root);

    std::vector<std::list<Vertex>> v_height; //element i should contain a list of all vertices in tree of height i
    int max_height = *max_element(std::begin(heights), std::end(heights));   
    for (int i = 0; i <= max_height; i++) {
        v_height.push_back({});
    }

    std::vector<std::list<Vertex>> v_colour; // element i should contain a list of all vertices in G with colour i
    int max_colour = boost::num_vertices(tree)-1; 

    for (int i = 0; i <= max_colour; i++) {
        v_colour.push_back({});
    }

    Graph::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(tree); vi != vi_end; ++vi) {
        Vertex v = *vi;
        v_height[heights[v]].push_back(v); //update height lookup for each vertex in tree
    }

    for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) { 
        Vertex v = *vi;
        v_colour[colour_G[v]].push_back(v); //update colour lookup for each vertex in G
    }
  
    // initialise count lookup to 0 for all vertices in G
    for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) {
        Vertex v = *vi;
        cpi_count[v] = 0;       
    }

    //iterate through each possible height for vertices in the tree
    for(int i = 0; i < v_height.size(); i++) {
        
        std::list<Vertex> all_height_i = v_height[i];
        // iterate through all vertices of height i
        for (Vertex v: all_height_i) {
            std::list<Vertex> candidates = v_colour[colour_tree[v]]; //v can only be mapped to vertices in G with the same colour
            for (Vertex w: candidates) {
                cpi_count[w] = 1;
                if (i != 0) {


                    //iterate through the children of v in tree
                    std::pair<DiGraph::out_edge_iterator, DiGraph::out_edge_iterator> outEdges = boost::out_edges(v, tree);
                    for (DiGraph::out_edge_iterator outEdge = outEdges.first; outEdge != outEdges.second; ++outEdge) {
                        Vertex c = boost::target(*outEdge, tree);
                        int current_child = 0;

                        //iterate through adjacent vertex to c in G
                        std::pair<Graph::adjacency_iterator, Graph::adjacency_iterator> neighbours = boost::adjacent_vertices(w, G);
                        for (Graph::adjacency_iterator neighbour = neighbours.first; neighbour != neighbours.second; ++neighbour) {
                            Vertex c_cand = *neighbour;
                            if (colour_G[c_cand] == colour_tree[c]) {
                                // c_cand can be mapped to by c
                                current_child += cpi_count[c_cand];
                            }                        
                        }


                        cpi_count[w] *= current_child;
                    }

                    if (i == v_height.size()-1) { //if on final iteration, i.e the root vertex
                        total += cpi_count[w]; //sum up all possible counts from the root vertex
                    }
                }

            }
        }   

        // after iteration i, all counts for vertices in G with the same colour as vertices in tree of height at most i are correct
    }

    return total;
}