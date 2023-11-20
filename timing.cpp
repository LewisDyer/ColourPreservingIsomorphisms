
# include <chrono>
# include <ctime>
# include "graph_utils.h"
# include "graph_generators.h"
# include "count_colourful_isomorphisms.h"
# include "performance_analysis.h"

using namespace boost;
using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Vertex =  graph_traits<DiGraph>::vertex_descriptor;
using Edge = graph_traits<DiGraph>::edge_descriptor;
using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
using Colours= std::map<Vertex, int>;
using ColourMap = associative_property_map<Colours>;

int main() {

    Graph GH = path(2);
    int root = 0;
    DiGraph H = makeRooted(GH, root);

    Graph G = star(50000);

    Colours col_H;

    // uniform colouring for H
    for(int i=0; i < boost::num_vertices(H); i++) {
        col_H[i] = i;
    }

    ColourMap colour_H(col_H);

    Colours col_G;

    for(int i=0; i < boost::num_vertices(G); i++) {
        if (i == 0) {
            col_G[i] = 0;
        } else {
            col_G[i] = 1; // central vertex has colour 0, all others have colour 1.
        }
    }

    // // uniform colouring for G
    // for(int i=0; i < boost::num_vertices(G); i++) {
    // col_G[i] = i % boost::num_vertices(H);
    // }

    ColourMap colour_G(col_G);

    auto now = std::chrono::system_clock::now();

    // Convert the time point to a time_t object
    std::time_t currentTime = std::chrono::system_clock::to_time_t(now);

    // Convert the time_t object to a tm structure
    std::tm* timeInfo = std::localtime(&currentTime);

    // Format the date and time as a string in ISO 8601 format
    std::ostringstream oss;
    oss << std::put_time(timeInfo, "%Y-%m-%dT%H-%M-%S");

    std::string timeString = oss.str();

    std::string filename = "outputs/" + timeString;


    time_tree_count(H, root, col_H, G, col_G, filename, "A path on 2 vertices", "A star on 50000 vertices with central vertex colour 0 and all other vertices colour 1", 100);

    // int k = 100; //no. of executions

    // for (int i=0; i < k; i++) {
    //     auto start = std::chrono::high_resolution_clock::now();
    //     int j = new_tree_count(H, 1, colour_H, G, colour_G);
    //     auto end = std::chrono::high_resolution_clock::now();

    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    //     std::cout << duration.count() << "ms \n";
    // }

    return 0;
}