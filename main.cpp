# include <iostream>
# include "count_colourful_isomorphisms.cpp"

# include "graph_utils.h"
# include "graph_generators.h"

int main() {

// uStar is rooted at vertex 1
// note root selection doesn't change the overall pattern, just how pattern matching is processed.
Graph uStar = star(4);
DiGraph H = makeRooted(uStar, 1);

Graph G = erdos_renyi(1000, 0.10);

Colours col_H;

// uniform colouring for H
for(int i=0; i < boost::num_vertices(H); i++) {
    col_H[i] = i;
}

ColourMap colour_H(col_H);

save_graph("H_coloured.dot", H);

Colours col_G;

// uniform colouring for G
for(int i=0; i < boost::num_vertices(G); i++) {
   col_G[i] = i % boost::num_vertices(H);
}

ColourMap colour_G(col_G);

save_graph("G_coloured.dot", G, colour_G);

int count = new_tree_count(H, 1, colour_H, G, colour_G);

std::cout << "NUMBER OF CPIs IS " << count << "\n";

}