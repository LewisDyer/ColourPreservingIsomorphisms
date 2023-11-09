
# include <iostream>
# include "count_colourful_isomorphisms.cpp"

# include "graph_utils.h"
# include "graph_generators.h"

int main() {

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