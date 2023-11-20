/***
 * Contains various timing methods for new_tree_count 
*/

# include <chrono>
# include <random>
# include "graph_utils.h"
# include "graph_generators.h"
# include "count_colourful_isomorphisms.h"

using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

/**
 * @param tree, the tree pattern being searched for
 * @param root the root vertex of the tree pattern
 * @param colour_H a property map defining a colourful k-vertex colouring of the tree
 * @param G an undirected data graph
 * @param colour_G a property map defining a k-vertex colouring of G (not necessarily colourful)
 * @param outputFile the filename to output results to
 * @param describe_tree A short user-provided description of the tree pattern.
 * @param describe_G A short user-provided description of the data graph
 * @param noRuns The number of times we run new_tree_count
 * 
 * Runs new_tree_count the specified number of times with the given parameters, outputting the results and summary statistics to outputFile.
*/
void time_tree_count(DiGraph tree, Vertex root, ColourMap colour_tree, Graph G, ColourMap colour_G, std::string outputFile, std::string describe_tree, std::string describe_G, int noRuns) {
    std::ofstream file;
    int total_runtime = 0;
    file.open(outputFile);
    file << describe_tree << "\n";
    file << describe_G << "\n";
    file << noRuns << " executions:\n";
    for (int i=0; i < noRuns; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        int j = new_tree_count(tree, root, colour_tree, G, colour_G);
        auto end = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);


        std::cout << duration.count() << "ms \n";
        file << duration.count() << "ms \n";
        total_runtime += duration.count();
    }

    long meanTime = total_runtime / noRuns;

    file << "Mean runtime: " << meanTime << "ms \n";

    file.close();
}

void time_with_random_colouring(DiGraph tree, Vertex root, ColourMap colour_tree, Graph G, std::string outputFile, std::string describe_tree, std::string describe_G, int noRuns) {
    std::ofstream file;
    int total_runtime = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> randColour(0, boost::num_vertices(tree)-1); // same number of colours as the tree
    Colours col_G;
    file.open(outputFile);
    file << describe_tree << "\n";
    file << describe_G << "\n";
    file << noRuns << " executions:\n";
    file << "with random colourings for G per execution\n";
    for (int i=0; i < noRuns; i++) {

        for(int i=0; i < boost::num_vertices(G); i++) {
            col_G[i] = randColour(gen);
        }
        
        ColourMap colour_G(col_G);
        auto start = std::chrono::high_resolution_clock::now();
        int j = new_tree_count(tree, root, colour_tree, G, colour_G);
        auto end = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);


        std::cout << duration.count() << "ms \n";
        file << duration.count() << "ms \n";
        total_runtime += duration.count();
    }

    long meanTime = total_runtime / noRuns;

    file << "Mean runtime: " << meanTime << "ms \n";

    file.close();
}