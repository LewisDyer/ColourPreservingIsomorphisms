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

// TREE-SPECIFIC IMPLEMENTATION
// r is current vertex being considered, mapped is the vertex that the root is mapped to in the current partial solution
// will need a wrapper for the root call
// Count(r, mapped):
//      if r is a leaf node:
//          return 1
//      else:
//          for each child of r, denoted c:
//              current_count = count(c, i) for each i adjacent to mapped with colour equal to r
//          return product of current_counts



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
# include "graph_generators.h"

using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
//using BagsMap = std::map<Vertex, std::set<Vertex>>;
//using DecompositionBags = associative_property_map<BagsMap>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

//# include "tree_decomposition.hpp"
//# include "nice_tree_decomposition.hpp"
//# include "graph_utils.h"
//# include "graph_generators.cpp"

/**
 * @param n, the number of vertices in the graph.
 * @return A path on n vertices.
*/
// template <class G>
// G path(int n) {

//     G path = DiGraph(n);

//     for(int i=0; i < n-1; i++) {
//         add_edge(i, i+1, path);
//     }

// 	return path;
// }

// template <class G>
// G star (int n) {
//     G star = Graph(n);

//     for (int i=1; i < n; i++) {
//         add_edge(0, i, star);
//     }

//     return star;
// }

/**
 * @param n, the number of vertices in the graph.
 * @return A clique on n vertices. 
*/
// template <class G>
// G clique(int n) {
//     G clique = Graph(n);
//     for(int i=0; i < n; i++) {
//         for (int j=0; j < i; j++) {
//             add_edge(i, j, clique);
//         }    
//     }

//     return clique;
// }

// namespace std {
// 	template <class T>
//     // Function to print out a set
// 	std::ostream& operator<< (std::ostream& os, const std::set<T> & in) {
// 		os << "[";
// 		for (auto it = in.begin(); it != in.end(); it++) {
//             if (it != in.begin()) os << ",";
//             os << *it;
//         }
//         os << "]";

// 		return os;
// 	}
// }

// void output_bags(Graph G, DecompositionBags bags_pm) {
//     for(int i=0; i < boost::num_vertices(G); i++) {
//         std::cout << i << ": " << bags_pm[i] << "\n";
//     }
// }

// Graph erdos_renyi(int n, double p) {
//     Graph g(n);

//     boost::random::mt19937 gen;
//     gen.seed(std::time(0));

//     boost::random::uniform_int_distribution<> dis(0, n-1);
//     boost::random::uniform_real_distribution<double> real_dis(0.0, 1.0);
//     for (int i = 0; i < n; ++i) {
//         for (int j = i + 1; j < std::min(i+1+(int)(n*p), n); ++j) {
//             if (real_dis(gen) < p) {
//                 boost::add_edge(i, j, g);
//             }
//         }
//     }

//     return g;
// }

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
// Vertex get_root(DiGraph tree) {
//     typedef typename graph_traits<DiGraph>::vertex_iterator vi;
//     for (std::pair<vi, vi> ver = boost::vertices(tree); ver.first != ver.second; ++ver.first) {
//         // Check the in-degree of each vertex
//         if (boost::in_degree(*ver.first, tree) == 0) {
//             // Return the first vertex with in-degree 0 found
//             return *ver.first;
//         }
//     }

//     return 0;
// }



// enum VertexType {leaf, join, forget, introduce};
// const char* node_type[4] = {"leaf", "join", "forget", "introduce"};

// /**
//  * @param v the vertex we're considering in the nice tree decomposition
//  * @param tree the nice tree decomposition containing v
//  * @param bags the property map of bags in the tree
//  * @return The type of vertex that v is in tree (either leaf, join, forget or introduce) 
// */
// VertexType getNiceType(const Vertex& v, const DiGraph& tree, DecompositionBags bags) {
//     // Count children
//     // If 0 children, return LEAF
//     // If 2 children, return JOIN
//     // Then compare bag of v and the bag of child of v
//     // If v's bag is bigger return FORGET, else return INTRODUCE

//     switch(boost::out_degree(v, tree)) {
//         case 0: { return leaf; break; }
//         case 2: { return join; break; }
//         case 1: {
//             typename graph_traits<DiGraph>::out_edge_iterator child, child_end;
//             boost::tie(child, child_end) = out_edges(v, tree);
//             VertexType type = (size(bags[v]) > size(bags[boost::target(*child, tree)])) ? introduce : forget;
//             return type;
//             break;
//         }
//         default: { return leaf; break; }
//     }
// }

// Graph induced_subgraph(Graph G, std::set<Vertex> V) {
//     Graph S(V.size());

//     for (Vertex v: V) {
//         Vertex v_sub = vertex(v, S);
//         auto adj_vertices = boost::adjacent_vertices(v, G);
//         for (auto v_target: boost::make_iterator_range(adj_vertices.first, adj_vertices.second)) {

//             if (std::find(V.begin(), V.end(), v_target) != V.end()) {
//                 if (!boost::edge(v_target,v,S).second) {add_edge(v, v_target, S); }        
//             }
            
//         }
//     }

//     return S;
// }

// struct iso_params {
//     std::map<Vertex, Vertex> H_index_map;
//     std::map<Vertex, Vertex> G_index_map;
//     ColourMap colour_H;
//     ColourMap colour_G;
// };

// /**
//  * @param v1 a vertex in H
//  * @param v2 a vertex in G
//  * @param H a subgraph of the pattern graph
//  * @param G a subgraph of the data graph
//  * @param params a struct containing H_index_map, sending the logical representatives of the subgraph of H to the graph representation, and G_index_map which does the opposite for G
//  * Also colour_H and colour_G, which map logical representatives of H and G respectively to their colours.
//  * @return true if v1 and v2 have the same colour, false if they don't.
// */
// bool colour_match(Vertex v1, Vertex v2, Graph H, Graph G, iso_params params) {

//     std::map<Vertex, Vertex> rev_H_index; // from graph rep of sub(H) to logical rep of sub(H);

//     for (auto pair: params.H_index_map) {
//             rev_H_index[pair.second] = pair.first;
//         }

//     int H_colour = params.colour_H[rev_H_index[v1]];
//     int G_colour = params.colour_H[params.G_index_map[v2]];

//     return (H_colour == G_colour);

// }

//  std::pair<Graph, std::map<Vertex, Vertex>> induced_subgraph(Graph G, std::set<Vertex> V, bool inverseMap) {

//     std::map<Vertex, Vertex> index_map;

//     if (V.empty()) { Graph H(0); return std::make_pair(H, index_map);}
//     Graph H;


//     int i = 0;
//     for (Vertex u:V) {

//         // std::cout << "add vertex " << u << "\n";
//         add_vertex(u, H);
//         index_map[u] = i;
//         ++i;
//     }

//     //for(auto v: index_map) {std::cout << v.first << " -> " << v.second << "\n";}
//     //std::cout << "----\n";

//     // std::cout << "# of vertices in graph: " << num_vertices(H) << "\n";

//     for (Vertex u: V) {
//         for (auto e: make_iterator_range(out_edges(u,G))) {
//             Vertex v = target(e, G);
//             //std::cout << "checking edge (" << index_map[u] << ", " << index_map[v] << ")\n";
//             //std::cout << "without index_map it's (" << u << ", " << v << ")\n";
            
//             if (find(V.begin(), V.end(), v) != V.end()) {
//                 if (u < v) {

//                 // std::cout << "add edge (" << index_map[u] << ", " << index_map[v] << ")\n";

//                 add_edge(index_map[u], index_map[v], H);

//             }
//             }
//         }
//     }

//     if (!inverseMap) {
//         return std::make_pair(H, index_map);

//     } else {
//         std::map<Vertex, Vertex> reverse_map;

//         for (auto pair: index_map) {
//             reverse_map[pair.second] = pair.first;
//         }

//         return std::make_pair(H, reverse_map);
//     }

    
//  }

// //     boost::graph_traits<Graph>::vertex_iterator vi, vi_end, next;
// //     if (boost::num_vertices(H) <= 1) {return H;}

// //     boost::tie(vi, vi_end) = boost::vertices(H);
// //     for(next=vi; vi != vi_end; vi=next) {
// //         if (boost::degree(*vi, H) == 0) {
// //             boost::remove_vertex(*vi, H);
// //         }
// //         ++next;
// //     }

// //     if (boost::degree(0, H) == 0) {
// //         boost::remove_vertex(0, H);
// //     }

// //     return H;
// // }

int count_tree(Vertex root, Vertex rootMapped, DiGraph tree, ColourMap colour_H, Graph G, ColourMap colour_G) {

    // Count(r, mapped):
    //      if r is a leaf node:
    //          return 1
    //      else:
    //          for each child of r, denoted c:
    //              current_count = count(c, i) for each i adjacent to mapped with colour equal to r
    //          return product of current_counts

    if (boost::out_degree(root, tree) == 0 ) { //leaf node
        return 1;
    }
    
    int combos = 1;
    // iterate through children of r

    typename graph_traits<DiGraph>::out_edge_iterator child, child_end;
    typename graph_traits<DiGraph>::out_edge_iterator target, target_end;

    for(boost::tie(child, child_end) = boost::out_edges(root, tree); child != child_end; ++child) {
        Vertex childVertex = boost::target(*child, tree);
        int childDegree = boost::degree(childVertex, tree);
        int childColour = colour_H[childVertex];
        int currentTotal = 0;
        //iterate through adjacent vertices of rootMapped, look for vertices with childColour
        
        std::pair<Graph::adjacency_iterator, Graph::adjacency_iterator> neighbors = boost::adjacent_vertices(rootMapped, G);

        // Iterate through adjacent vertices
        for (Graph::adjacency_iterator neighbor = neighbors.first; neighbor != neighbors.second; ++neighbor) {
            Graph::vertex_descriptor adjacent_vertex = *neighbor;
            int targetColour = colour_G[adjacent_vertex];
            if (targetColour == childColour) {
                int targetDegree = boost::degree(adjacent_vertex, G);
                if (targetDegree >= childDegree) { // degree filtering to avoid checking impossible matches
                currentTotal += count_tree(childVertex, adjacent_vertex, tree, colour_H, G, colour_G);
                } else {
                    std::cout << "avoided due to degree filtering \n";
                }
            }
        }

        if (currentTotal == 0) {return 0;}

        combos *= currentTotal;

    }
    
    return combos;

}

/**
 * @param tree the pattern graph being considered (must be a tree on k vertices)
 * @param root The root of the tree
 * @param colour_H the colourful vertex k-colouring of the vertices of the tree
 * @param G the data graph
 * @param colour_G an arbitrary k-colouring of vertices of G
 * @return The number of colour-preserving isomorphisms between H and subgraphs of G.
*/
int tree_count(DiGraph tree, Vertex root, ColourMap colour_H, Graph G, ColourMap colour_G) {

    //Vertex root = get_root(tree);
    int targetColour = colour_H[root];
    int rootDegree = boost::degree(root, tree);

    // search for targetColour in G to consider starting points

     int total = 0;
            typedef typename graph_traits<Graph>::vertex_iterator iter_v;
            for (std::pair<iter_v, iter_v> p = vertices(G); p.first != p.second; ++p.first) {
                if (colour_G[*p.first] == targetColour && *p.first != boost::num_vertices(G)) {
                    if (boost::degree(*p.first, G) >= rootDegree) {
                    total += count_tree(root, *p.first, tree, colour_H, G, colour_G);
                }
                }
            }


    return total;
}

// Function to calculate vertex heights in a tree
void calculateHeights(const DiGraph& tree, Vertex v, std::vector<int>& heights) {
    //std::cout << "CALLED FOR VERTEX " << v << "\n";
    int height = 0; // Initialize the height for the current vertex
    
    // Loop through all adjacent vertices
    for (auto it = boost::adjacent_vertices(v, tree); it.first != it.second; ++it.first) {
        int adjacent_vertex = *it.first;
        calculateHeights(tree, adjacent_vertex, heights);
        height = std::max(height, heights[adjacent_vertex] + 1); // Update flavor
    }
    
    //std::cout << "VERTEX " << v << " HAS HEIGHT " << height << "\n";
    heights[v] = height; // Set the flavor of the current vertex
}

std::vector<int> getAllHeights(const DiGraph& tree, Vertex root) {
    std::vector<int> heights(boost::num_vertices(tree), 0);
    calculateHeights(tree, root, heights);
    return heights;
}

int new_tree_count(DiGraph tree, Vertex root, ColourMap colour_H, Graph G, ColourMap colour_G) {
    // std::cout<<"Starting the algorithm inner\n";
    int total = 0;
    Colours cpi_count; //stores counts of partial solutions for all vertices in G
    typedef typename graph_traits<Graph>::vertex_iterator iter_v;
    std::vector<int> heights = getAllHeights(tree, root);
    // std::cout << "HEIGHTS: [";
    // for (int h: heights) {
    //     std::cout << h << " ";
    // }
    // std::cout << "]\n";



    std::vector<std::list<Vertex>> v_height; // element i contains a list of vertices in T of height i
    int max_height = *max_element(std::begin(heights), std::end(heights));
    
    for (int i = 0; i <= max_height; i++) {
        v_height.push_back({});
    }

    std::vector<std::list<Vertex>> v_colour; // element i contains a list of vertices with colour i
    int max_colour = boost::num_vertices(tree)-1;

    for (int i = 0; i <= max_colour; i++) {
        v_colour.push_back({});
    }

    // std::cout << v_height.size();
    // std::cout << v_colour.size();

    // std::cout<<"v_colour map\n";

    // for (std::list<Vertex> colour: v_colour) {



    //     std::cout<<"[";
    //     for (Vertex v: colour) {
    //         std::cout << v << " ";
    //     }
    //     std::cout<<"]\n";
    // }


    //std::cout<<"Starting vertex traversal\n";
    Graph::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(tree); vi != vi_end; ++vi) { //think this traversal might be odd
        Vertex v = *vi;
        //std::cout<<"Vertex has height " << heights[v] << "\n";
        //std::cout<<"And colour " << colour_H[v] << "\n";
        v_height[heights[v]].push_back(v);
    }

    for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) { //think this traversal might be odd
        Vertex v = *vi;
        //std::cout << "Vertex " << v << " has colour " << colour_G[v] << "\n";
        v_colour[colour_G[v]].push_back(v);
    }



    // std::cout<<"v_height map\n";

    // for (std::list<Vertex> height: v_height) {

    //     std::cout<<"[";
    //     for (Vertex v: height) {
    //         std::cout << v << " ";
    //     }
    //     std::cout<<"]\n";
    // }

    // std::cout<<"v_colour map\n";

    // for (std::list<Vertex> colour: v_colour) {

    //     std::cout<<"COLOUR SIZE:\n";

    //     std::cout << colour.size() << "\n";

    //     std::cout<<"[";
    //     for (Vertex v: colour) {
    //         std::cout << v << " ";
    //     }
    //     std::cout<<"]\n";
    // }

    // std::cout << boost::num_vertices(G);

    
    for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) {
        Vertex v = *vi;
        //std::cout << "INITIALISING FOR " << v << "\n";
        cpi_count[v] = 0;
        
    }
        for(int i = 0; i < v_height.size(); i++) {
            std::list<Vertex> all_height_i = v_height[i];
            for (Vertex v: all_height_i) {
                std::list<Vertex> candidates = v_colour[colour_H[v]];
                for (Vertex w: candidates) {
                    //std::cout << "CHECKING CANDIDATE " << w << "\n";
                    cpi_count[w] = 1;
                    if (i != 0) {

                        std::pair<DiGraph::out_edge_iterator, DiGraph::out_edge_iterator> outEdges = boost::out_edges(v, tree);

                        for (DiGraph::out_edge_iterator outEdge = outEdges.first; outEdge != outEdges.second; ++outEdge) {
                            Vertex c = boost::target(*outEdge, tree);
                            int current_child = 0;
                            //std::cout << "RESET CURRENT CHILD\n";

                            std::pair<Graph::adjacency_iterator, Graph::adjacency_iterator> neighbours = boost::adjacent_vertices(w, G);
                            for (Graph::adjacency_iterator neighbour = neighbours.first; neighbour != neighbours.second; ++neighbour) {
                                Vertex c_cand = *neighbour;
                                if (colour_G[c_cand] == colour_H[c]) {
                                    //std::cout << "C' IS VERTEX " << c_cand << "\n";
                                    //std::cout << "INCREMENT CURRENT CHILD BY " << cpi_count[c_cand] << "\n";
                                    current_child += cpi_count[c_cand];
                                }                        
                            }

                            //std::cout << "MULTIPLY BY " << current_child << "\n";

                            cpi_count[w] *= current_child;
                        }

                        if (i == v_height.size()-1) {
                            total += cpi_count[w];
                        }
                    }

                    //std::cout << "CANDIDATE " << w << " HAS CPI COUNT " << cpi_count[w] << "\n";
                }
            }   
        }

    return total;
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
 * @param isCounting True for counting, false for just deciding.
 * @return The number of colour-preserving isomorphisms between subgraphs of G and H_y that extend colour-preserving isomorphisms between K and X_y. If no such isomorphisms exist between K and X_y, return 0.
*/
// int colourful_count(Vertex root, DiGraph tree, std::set<Vertex> K, Graph H, Graph G, DecompositionBags bags,  ColourMap colour_H, ColourMap colour_G, bool isCounting) {
//     // Get the vertices in the root bag
//     // Get the subgraph of H induced by these vertices

//     //std::cout <<"=====NEW VERTEX=====\n";

//     std::set<Vertex> root_vertices = bags[root];

//    //std::cout << "root is " << root_vertices << "\n";

//     //std::cout << "K is " << K << "\n";

//     std::map<Vertex, Vertex> index_map_H, index_map_G;
//     Graph Xy, Gk;


//     boost::tie(Xy, index_map_H) = induced_subgraph(H, root_vertices, false);
//     boost::tie(Gk, index_map_G) = induced_subgraph(G, K, true);

    

//     //save_graph("Test_Xy.dot", Xy);

//     //save_graph("Test_Gk.dot", Gk);

//     std::vector<Vertex> iso(num_vertices(H));

//     //std::cout << root_vertices << " in H iso. to " << K << " in G? ";

//     bool areIsomorphic;

//     // trying to do isomorphism checks with empty graphs causes early termination
//     if ((num_vertices(Xy) == 0) && (num_vertices(Gk) == 0)) {
//         areIsomorphic = true;
//     } else if ((num_vertices(Xy) == 0) || (num_vertices(Gk) == 0)) {
//         areIsomorphic = false;
//     } else if ((num_vertices(Xy) == 1) && (num_vertices(Gk) == 1)) {
//         // NOTE will need to refactor this for colour checking too
//         areIsomorphic = true;
//         Vertex y = *(root_vertices.begin());
//         Vertex k = *(K.begin());
//         iso[y] = k;
//     } else {
//         // NOTE this doesn't check colour-preserving isomorphisms yet!
//         //areIsomorphic = boost::isomorphism(Xy, Gk, isomorphism_map(make_iterator_property_map(iso.begin(), get(vertex_index, Xy))));
//         iso_params params;
//         params.G_index_map = index_map_G;
//         params.H_index_map = index_map_H;
//         params.colour_G = colour_G;
//         params.colour_H = colour_H;

//         areIsomorphic = boost::isomorphism(Xy, Gk, boost::vertex_invariant(boost::bind(colour_match, _1, _2, boost::ref(Xy), boost::ref(Gk), boost::ref(params))));
//     }

//     if (!areIsomorphic) { return 0; }

//    //std::cout << "Yes\n=====\n";

//     //for (int i=0; i < iso.size(); i++) { std::cout << i << " maps to " << iso[i] << "\n"; }
    
    

//     //for (int i=0; i < iso.size(); i++) { std::cout << i << " maps to " << index_map_G[iso[i]] << "\n"; }
     
//     // // Check the type of the root node in the tree ("leaf"/"join"/"forget"/"introduce")

//      VertexType root_type = getNiceType(root, tree, bags);

//      //std::cout << "root type is " << node_type[root_type] << "\n" ;

//      switch(root_type) {
//         case leaf: {
//             return 1;
//             break;
//         }
//         case join: {
            
//             int count = 0;
//             Vertex child1, child2;

//             typename graph_traits<DiGraph>::out_edge_iterator child, child_end;

//             boost::tie(child, child_end) = boost::out_edges(root,tree);
//             child1 = boost::target(*child, tree);

//             ++child;
//             child2 = boost::target(*child, tree);

//             if (isCounting) {
//             return (colourful_count(child1, tree, K, H, G, bags, colour_H, colour_G, isCounting) * colourful_count(child2, tree, K, H, G, bags, colour_H, colour_G, isCounting));
//             } else {
//                 if (colourful_count(child1, tree, K, H, G, bags, colour_H, colour_G, isCounting) > 0) {
//                     return colourful_count(child2, tree, K, H, G, bags, colour_H, colour_G, isCounting);
//                 } else {
//                     return 0;
//                 }
//             }
//             break;
//         }
//         case forget: {

//             typename graph_traits<DiGraph>::out_edge_iterator child, child_end;
//             boost::tie(child, child_end) = boost::out_edges(root, tree);

//             Vertex child1 = boost::target(*child, tree);

//             //Get the colour of the new vertex in this node
//             std::set<Vertex> diff;
//             std::set_difference(get(bags, child1).begin(), get(bags, child1).end(), get(bags, root).begin(), get(bags, root).end(), std::inserter(diff, diff.begin()));
//             Vertex new_v = *(diff.begin());
//             int colour = colour_H[new_v];

//             // Iterate over vertices in G and check their colour - branch on all vertices with the same colour as the new vertex.

//             int total = 0;
//             typedef typename graph_traits<Graph>::vertex_iterator iter_v;
//             for (std::pair<iter_v, iter_v> p = vertices(G); p.first != p.second; ++p.first) {
//                 if (colour_G[*p.first] == colour && *p.first != boost::num_vertices(G)) {
//                     //std::cout << "try mapping " << new_v << " to " << *p.first << "\n";
//                     std::set<Vertex> new_K;
//                     std::copy(K.begin(), K.end(), std::inserter(new_K, new_K.begin()));
//                     new_K.insert(*p.first);
//                     if (isCounting) { total += colourful_count(child1, tree, new_K, H, G, bags, colour_H, colour_G, isCounting); }
//                     else if (colourful_count(child1, tree, new_K, H, G, bags, colour_H, colour_G, isCounting) > 0) {return 1;}
//                 }
//             }

//             return total;
//             break;
//         }
//         case introduce: {

//             typename graph_traits<DiGraph>::out_edge_iterator child, child_end;
//             boost::tie(child, child_end) = boost::out_edges(root, tree);

//             Vertex child1 = boost::target(*child, tree);

//             std::set<Vertex> diff;
//             std::set_difference(get(bags, root).begin(), get(bags, root).end(), get(bags, child1).begin(), get(bags, child1).end(), std::inserter(diff, diff.begin()));
//             Vertex forgotten_v = *(diff.begin());


//             // Given a vertex in Xy: get its representation in the subgraph variable, then the vertex it's mapped to in the representation of GK,
//             // then the actual logical representation of GK.
//             int prev_target = index_map_G[iso[index_map_H[forgotten_v]]];
//             // Not the most efficient, but size of this is bounded by size of H so not too bad
//             // for(int i=0; i < iso.size(); i++) {
//             //     if(iso[i] == forgotten_v) {
//             //         prev_target = iso[i];
//             //         break;
//             //     }
//             // }

//             //std::cout << "dropping vertex " << forgotten_v << " which was mapped to " << prev_target << "\n";
//             std::set<Vertex> new_K;
//             std::copy(K.begin(), K.end(), std::inserter(new_K, new_K.begin()));
//             new_K.erase(prev_target);
//             return colourful_count(child1, tree, new_K, H, G, bags, colour_H, colour_G, isCounting);
//             break;
//         }
//         default: {
//             std::cout << "default type (you should never see this!)\n";
//             return 0;
//             break;
//         }
//     }

//     // If "leaf", return 1
//     // If "join", return product of counts rooted children, use the same K
//     // If "forget", get the colour of the new node added in the child. Then iterate
//     // over all nodes in G with the same colour, and branch here to recurse rooting on the child with all possible extensions to K.
//     // If "introduce", find the vertex that was removed in the child, find its colour, then remove the vertex in K with that colour and recurse from the child.


//     return 1;
// }

/**
 * @param H the pattern graph.
 * @param G the data graph.
 * @param colour_H the property map containing the vertex colours of H, is a colourful k-colouring
 * @param colour_G the property map containing the vertex colours of G, is an arbitrary k-colouring
 * @param isCounting True if we want to count CPIs, False if we just want to check if any exist.
 * @return The number of subgraphs of G isomorphic to H and preserving vertex colourings
*/
// int count_colour_preserving_isomorphisms(Graph H, Graph G, ColourMap colour_H, ColourMap colour_G, bool isCounting) {
//     // get a tree decomposition of H
//     // make it nice, and get its root
//     // inner recursive function (see colourful_count)

//     DiGraph dec;
//     BagsMap bags_m;
//     DecompositionBags bags_pm(bags_m);
    
//     bool x = boost::tree_decomposition(H, dec, bags_pm, 1);

//     //std::cout << "initial decomposition fine\n";

//     DiGraph tree;
//     BagsMap nice_bags_m;
//     DecompositionBags bags(nice_bags_m);

//     Vertex root = boost::nice_tree_decomposition(dec, bags_pm, tree, bags);

//     //("test_nice_dec.dot", tree, bags);

//     // for(int i=0; i < boost::num_vertices(tree); i++) {
//     //     std::cout << "vertex " << i << " has node type " << node_type[getNiceType(i, tree, bags)] << "\n";
//     // }

//     //root always has empty bags, so K starts empty
//     std::set<Vertex> K;
//     return colourful_count(root, tree, K, H, G, bags, colour_H, colour_G, isCounting);
//     //return 1;
// }

void test() {

// Graph h = path(5);
// std::set<Vertex> V;
// V.insert(3);
// V.insert(4);

// std::cout<<"is this working?\n";

// Graph g = induced_subgraph(h,V);
// std::cout<<"try saving\n";
// save_graph("WHERE_IS_G_INDUCED.dot", g);
// save_graph("test_H.dot", h);

std::cout << "start method\n";

std::cout <<"with new algorithm\n";

Graph uStar = star(4);
DiGraph H = makeRooted(uStar, 1);

//DiGraph H = path<DiGraph>(4);

//Graph G = uPath<Graph>(10);
Graph G = erdos_renyi(1000, 0.10);

//std::cout << "made random graph\n";
//Graph G = path<Graph>(10);

//save_graph("random_G.dot", G);



Colours col_H;

for(int i=0; i < boost::num_vertices(H); i++) {
    col_H[i] = i;
}

ColourMap colour_H(col_H);

save_graph("H_coloured.dot", H);

Colours col_G;

boost::random::mt19937 gen;
gen.seed(std::time(0));

boost::random::uniform_int_distribution<> dis(0, boost::num_vertices(H)-1);

for(int i=0; i < boost::num_vertices(G); i++) {
   //col_G[i] = dis(gen);
   col_G[i] = i % boost::num_vertices(H);
}

// for(int i=0; i < boost::num_vertices(G); i++) {
//     if (i <=4) {
//         col_G[i] = i;
//     } else{
//         col_G[i] = 5 - std::abs(i-4);
//     }
   
// }

// std::set<Vertex> set_h, set_g;
// set_h.insert(0);
// set_h.insert(1);

// set_g.insert(8);
// set_g.insert(9);

// std::map<Vertex, Vertex> index_map;
// Graph h_ind, g_ind;

// //boost::tie(h_ind, index_map) = induced_subgraph(H, set_h);
// boost::tie(g_ind, index_map) = induced_subgraph(G, set_g);

// //save_graph("h_ind.dot", h_ind);
// save_graph("g_ind.dot", g_ind);

ColourMap colour_G(col_G);

save_graph("G_coloured.dot", G, colour_G);

//save_graph("H_colour.dot", H, colour_H);

//save_graph("G_colour.dot", G, colour_G);



// int count = count_colour_preserving_isomorphisms(H, G, colour_H, colour_G, false);

// int count = tree_count(H, 1, colour_H, G, colour_G);

// std::cout << "NUMBER OF CPIs IS " << count << "\n";

std::cout << "start to run new alg\n";

int count = new_tree_count(H, 1, colour_H, G, colour_G);

std::cout << "new alg terminating\n";

std::cout << "NUMBER OF NEW CPIs IS " << count << "\n";

}




// int main() { test();}