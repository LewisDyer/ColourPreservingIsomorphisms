#ifndef BOOST_GRAPH_TREE_DECOMPOSITION_HPP
#define BOOST_GRAPH_TREE_DECOMPOSITION_HPP

#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/tuple/tuple.hpp>

/**
 * @file tree_decomposition.hpp
 *
 * The \p tree_decomposition() constructs the tree decomposition of graph \p g, that will be stored in the graph \p d. Bags for each node will be stored in property map \p bags. Decomposition itself will be of width at most 3k+4 - if such decomposition is not possible to construct, the function return value will be false, otherwise true.
 */

namespace boost {
    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    namespace detail {
        #ifndef BOOST_GRAPH_INSERT_COLLECTION
        #define BOOST_GRAPH_INSERT_COLLECTION
        /**
         * Concept check of the concept InsertCollection
         * @tparam C Type of container that is being checked
         */
        template <class C>
        struct InsertCollection : Collection<C> {
        public:
            using Container_value_type = typename C::value_type;

            BOOST_CONCEPT_USAGE(InsertCollection)
            {
                x.insert(x.begin(), e);
            }
        private:
            C x;
            Container_value_type e;
        };
        #endif // BOOST_GRAPH_INSERT_COLLECTION

        /**
         * Retrieves vector containing components of graph g
         * @tparam Graph    Type of the graph
         * @param g         Graph from which all of its components will be obtained
         * @return Vector containing pointers to components of g
         */
        template <class Graph>
        std::vector<subgraph<Graph>*> get_components(subgraph<Graph> & g) {
            //vertex descriptor that is used to address vertices in graph
            using Graph_vertex = typename graph_traits<subgraph<Graph>>::vertex_descriptor;
            //type of value, that represents number of vertices in graph
            using Graph_vertex_size_type = typename graph_traits<subgraph<Graph>>::vertices_size_type;

            std::map<Graph_vertex, Graph_vertex_size_type> c;
            Graph_vertex_size_type num = connected_components(g, make_assoc_property_map(c));
            std::vector<subgraph<Graph>*> components(num);
            for (Graph_vertex_size_type i = 0; i < num; i++) components[i] = & g.create_subgraph();
            for (const Graph_vertex v : make_iterator_range(vertices(g))) add_vertex(g.local_to_global(v), *components[c[v]]);

            return components;
        }

        /**
         * Splits set s into two sets based on c_number
         * @tparam T        Type of elements of the set s
         * @param s         Set of elements that will be split
         * @param c_number  Seed number that will decide how the set s will be separated
         * @return tuple of two non-empty sets, that were created by splitting set \p s
         */
        template<class T>
        tuple<std::set<T>, std::set<T>> split_set(std::set<T> s, unsigned long long c_number) {
            //std::set of graph descriptors
            using Graph_vertex_set = std::set<T>;

            Graph_vertex_set a, b;
            unsigned long long set_size = s.size();
            for (unsigned long long i = 0; i < set_size; i++) {
                auto it = s.begin();
                if (c_number % 2 == 0) a.insert(*it);
                else b.insert(*it);
                s.erase(it);
                if (c_number != 0) c_number = c_number >> 1;
            }

            return make_tuple<Graph_vertex_set, Graph_vertex_set>(a, b);
        }

        /**
         * Finds vertex separation of graph parent using the max flow algorithm
         * @tparam Graph    Type of the graph
         * @param split     Tuple of two sets, where one will be source and the other sink
         * @param parent    Graph on which will be applied max flow algorithm
         * @param k         Parameter that limits size of desired vertex separation
         * @return Set of vertices that represents vertex separation of graph parent or empty set, if such separation doesn't exist
         */
        template <class Graph>
        std::set<typename graph_traits<Graph>::vertex_descriptor>
        max_flow_sep(const tuple<std::set<typename graph_traits<Graph>::vertex_descriptor>, std::set<typename graph_traits<Graph>::vertex_descriptor>> & split,
                     subgraph<Graph> & parent,
                     unsigned long k) {
            //vertex descriptor that is used to address vertices in graph
            using Graph_vertex = typename graph_traits<Graph>::vertex_descriptor;
            //edge descriptor that is used to address edges in graph
            using Graph_edge = typename graph_traits<Graph>::edge_descriptor;
            //std::set of graph descriptors
            using Graph_vertex_set = std::set<Graph_vertex>;
            using DGTraits = adjacency_list_traits<vecS, vecS, directedS>;
            //type of helper graph, on which will be max flow algorithm applied
            using DGraph = adjacency_list<listS,
                                          vecS,
                                          directedS,
                                          property<vertex_name_t, std::string>,
                                          property<edge_capacity_t, long, property<edge_residual_capacity_t, long, property<edge_reverse_t, DGTraits::edge_descriptor>>>>;
            //vertex descriptor that is used to address vertices in helper graph
            using DGraph_vertex = DGTraits::vertex_descriptor;
            //edge descriptor that is used to address edges in helper graph
            using DGraph_edge = DGTraits::edge_descriptor;

            Graph_vertex_set flow_sep;
            DGraph helper_graph;
            property_map<DGraph, edge_capacity_t>::type capacity_map = get(edge_capacity, helper_graph);
            property_map<DGraph, edge_reverse_t>::type rev_map = get(edge_reverse, helper_graph);
            std::map<Graph_vertex, std::pair<DGraph_vertex, DGraph_vertex>> mapper;

            //converting vertices to double vertices
            for (const Graph_vertex v : make_iterator_range(vertices(parent))) {
                DGraph_vertex v1 = add_vertex(helper_graph);
                DGraph_vertex v2 = add_vertex(helper_graph);
                mapper[v] = std::make_pair(v1, v2);
            }
            //recreating edges (now directed with capacity set)
            for (const Graph_edge e : make_iterator_range(edges(parent))) {
                DGraph_edge e_new = add_edge(mapper[source(e, parent)].second, mapper[target(e, parent)].first, helper_graph).first;
                DGraph_edge e_new_rew = add_edge(mapper[target(e, parent)].second, mapper[source(e, parent)].first, helper_graph).first;
                capacity_map[e_new] = 2;
                capacity_map[e_new_rew] = 2;
                rev_map[e_new] = e_new_rew;
                rev_map[e_new_rew] = e_new;
            }
            //creating vertex-edges
            for (const auto & curr : mapper) {
                DGraph_edge e = add_edge(curr.second.first, curr.second.second, helper_graph).first;
                DGraph_edge er = add_edge(curr.second.second, curr.second.first, helper_graph).first;
                capacity_map[e] = 1;
                capacity_map[er] = 0;
                rev_map[e] = er;
                rev_map[er] = e;
            }
            //creating source of flow
            DGraph_vertex src = add_vertex(helper_graph);
            for (const Graph_vertex curr : get<0>(split)) {
                DGraph_edge e = add_edge(src, mapper[curr].first, helper_graph).first;
                DGraph_edge er = add_edge(mapper[curr].second, src, helper_graph).first;
                capacity_map[e] = 2;
                capacity_map[er] = 0;
                rev_map[e] = er;
                rev_map[er] = e;
            }
            //creating sink of flow
            DGraph_vertex dest = add_vertex(helper_graph);
            for (const Graph_vertex curr : get<1>(split)) {
                DGraph_edge e = add_edge(mapper[curr].second, dest, helper_graph).first;
                DGraph_edge er = add_edge(dest, mapper[curr].first, helper_graph).first;
                capacity_map[e] = 2;
                capacity_map[er] = 0;
                rev_map[e] = er;
                rev_map[er] = e;
            }

            property_map<DGraph, edge_residual_capacity_t>::type res_capacity_map = get(edge_residual_capacity, helper_graph);
            std::vector<default_color_type> color(num_vertices(helper_graph));
            std::vector<DGTraits::edge_descriptor> pred(num_vertices(helper_graph));

            unsigned long flow = edmonds_karp_max_flow(helper_graph, src, dest, capacity_map, res_capacity_map, rev_map, &color[0], &pred[0]);
            if (flow > k + 1) return flow_sep;
            //iterate over mapper and insert white vertices into result set
            for (const auto & curr : mapper) {
                if (color[curr.second.first] != color[curr.second.second]) flow_sep.insert(curr.first);
            }

            return flow_sep;
        }

        /**
         * Finds vertex separation of set s such that size of it doesn't exceed k+1
         * @tparam Graph    Type of the graph
         * @param s         Set of vertices from parent
         * @param parent    Graph on which  max flow algorithm will be applied
         * @param k         Parameter that limits size of desired vertex separation
         * @return Set of vertices that represents vertex separation of set s (with size at most k+1) or empty set, if such separation doesn't exist
         */
        template <class Graph>
        std::set<typename graph_traits<Graph>::vertex_descriptor>
        get_separation(const std::set<typename graph_traits<Graph>::vertex_descriptor> & s,
                       subgraph<Graph> & parent,
                       unsigned long k) {
            //vertex descriptor that is used to address vertices in graph
            using Graph_vertex = typename graph_traits<Graph>::vertex_descriptor;
            //std::set of graph descriptors
            using Graph_vertex_set = std::set<Graph_vertex>;

            //each split must have at least k+2 elements
            if (s.size() < (2*k +4)) return Graph_vertex_set();
            //max_c_num is 2^(s.size - 1)
            unsigned long long max_c_num = 1 << (s.size() - 1);
            //generate one by one all combinations of separation of set s into two sets with at least k+2 elements
            for (unsigned long long i = 1; i < max_c_num; i++) {
                auto split = split_set(s, i);
                if (get<0>(split).size() < k + 2
                    || get<1>(split).size() < k + 2) continue;
                Graph_vertex_set separation = max_flow_sep(split, parent, k);
                if (!separation.empty()) return separation;
            }

            return Graph_vertex_set();
        }

        /**
         * Retrieves set of vertices, that are neighbours of vertices of graph target inside its parent graph
         * @tparam Graph            Type of the graph
         * @param target            Graph containing vertices, whose neighbours will be returned
         * @param including_target  True, if vertices of target should be included, false if not
         * @return Set of vertices that are neighbours to vertices of target
         */
        template <class Graph>
        std::set<typename graph_traits<Graph>::vertex_descriptor>
        get_neighbour_vertices(subgraph<Graph> & target, bool including_target) {
            //vertex descriptor that is used to address vertices in graph
            using Graph_vertex = typename graph_traits<Graph>::vertex_descriptor;
            //std::set of graph descriptors
            using Graph_vertex_set = std::set<Graph_vertex>;

            Graph_vertex_set res;
            //iterate over every vertex of target
            for (const Graph_vertex v : make_iterator_range(vertices(target))) {
                //insert vertex from target itself (we may delete it later on)
                res.insert(target.local_to_global(v));
                //get all neighbors and insert them into set
                auto neighbours = adjacent_vertices(target.local_to_global(v), target.root());
                for (const Graph_vertex neighbour : make_iterator_range(neighbours)) res.insert(neighbour);
            }
            //if we don't want to include target vertices, we delete them
            if (!including_target) {
                for (const Graph_vertex v : make_iterator_range(vertices(target))) res.erase(target.local_to_global(v));
            }

            return res;
        }

        /**
         * Constructs part of tree decomposition of graph \p parent, where following conditions are met:
         * 1. s ⊂ s1 ⊆ w
         * 2. |s1| ≤ 4k + 5
         * 3. every connected component of parent[w \ s1] is adjacent to at most 3k+4 vertices of V(s1)
         * where s1 is root bag of current call of decomposition()
         * @tparam Graph            Type of the graph
         * @tparam Decomposition    Type of the decomposition of Graph
         * @tparam Bags             Type of the bags of Decomposition
         * @param parent            Graph whose tree decomposition will be constructed
         * @param decomposition     Graph in which will be the decomposition stored
         * @param bags              Property map containing bags for each node of d
         * @param w                 Set of vertices of parent
         * @param s                 Set of vertices of parent
         * @param k                 Parameter, which will limit width of decomposition bags
         * @return Tuple, where first type is vertex descriptor of root of current calls decomposition tree and the second type is if decomposition was found (true) or not (false)
         */
        template <class Graph, class Decomposition, class Bags>
        tuple<typename graph_traits<Decomposition>::vertex_descriptor, bool>
        decompose(subgraph<Graph> & parent,
                  Decomposition & decomposition,
                  Bags & bags,
                  const std::set<typename graph_traits<Graph>::vertex_descriptor> & w,
                  const std::set<typename graph_traits<Graph>::vertex_descriptor> & s,
                  unsigned long k) {
            //vertex descriptor that is used to address vertices in graph
            using Graph_vertex = typename graph_traits<Graph>::vertex_descriptor;
            //vertex descriptor that is used to address vertices in decomposition graph
            using Decomposition_vertex = typename graph_traits<Decomposition>::vertex_descriptor;
            //std::set of graph descriptors
            using Graph_vertex_set = std::set<Graph_vertex>;

            Graph_vertex_set s1 = s;
            if (s.size() < (3*k + 4)) {
                //find vertex u, such that u is in W \ S
                for (const Graph_vertex u : w) {
                    if (s1.size() == (3*k + 4))  break;
                    s1.insert(u);
                }
            }
            else {
                subgraph<Graph> & w_graph = parent.create_subgraph(w.begin(), w.end());
                Graph_vertex_set separation = detail::get_separation<Graph>(s, w_graph, k);
                if (separation.empty()) return make_tuple(NULL, false);
                for (const Graph_vertex u : separation) s1.insert(u);
            }

            //constructing root bag of current call of decompose()
            Decomposition_vertex root_bag = add_vertex(decomposition);
            for (const Graph_vertex u : s1) bags[root_bag].insert(bags[root_bag].begin(), u);
            //constructing W \ S1 and its graph
            Graph_vertex_set d;
            set_difference(w, s1, std::inserter(d, d.begin()));
            subgraph<Graph> & d_graph = parent.create_subgraph(d.begin(), d.end());
            //get all components of d_graph
            std::vector<subgraph<Graph>*> d_components = detail::get_components(d_graph);
            //for every component we will recursively call decompose
            for (auto component : d_components) {
                tuple<Decomposition_vertex, bool> res = detail::decompose(
                        parent,
                        decomposition,
                        bags,
                        detail::get_neighbour_vertices(*component, true),
                        detail::get_neighbour_vertices(*component, false),
                        k);
                if (!get<1>(res)) return res;
                //connect bag with current root
                add_edge(get<0>(res), root_bag, decomposition);
            }

            return make_tuple(root_bag, true);
        }
    } // namespace detail
    #endif //DOXYGEN_SHOULD_SKIP_THIS

    /**
     * @tparam     Graph            Type of the graph.
     * @tparam     Decomposition    Type of the decomposition.
     * @tparam     Bags             Type of the property map representing bags of \p Decomposition.
     * @param[in]  g                An undirected graph. The graph type must be a model of VertexListGraph.
     * @param[out] d                An undirected graph in which will be stored the tree decomposition of \p g.
     * @param[out] bags             The bags property map. The type must be a model of a mutable Readable Property Map. The key type of the map must be a vertex descriptor of \p d. The value type of the map must be a model of Mutable_Container. It also must implement function insert([position/hint iterator], [element]). Container value type must be vertex descriptor of \p g.
     * @param[in]  k                Parameter, which defines the maximum width (3k+4) of the tree decomposition.
     * @return True if tree decomposition of graph g exists (with width at most 3k+4), false if not.
     */
    template <class Graph, class Decomposition, class Bags>
    bool tree_decomposition(Graph & g, Decomposition & d, Bags & bags, unsigned long k) {
        //vertex descriptor that is used to address vertices in graph
        using Graph_vertex = typename graph_traits<Graph>::vertex_descriptor;
        //vertex descriptor that is used to address vertices in decomposition graph
        using Decomposition_vertex = typename graph_traits<Decomposition>::vertex_descriptor;
        //type of value (container) of Property Map
        using Bag_value_type = typename property_traits<Bags>::value_type;
        //type of container value, that is stored as value in Property Map
        using Bag_inside_type = typename Bag_value_type::value_type;

        //Graph must be model of Vertex List Graph
        BOOST_CONCEPT_ASSERT((VertexListGraphConcept<Graph>));
        //each bag of Bags container must be a model of Mutable_Container
        BOOST_CONCEPT_ASSERT((detail::InsertCollection<Bag_value_type>));
        //Bags must be ReadWrite Property Map (with Decomposition_vertex as a key)
        BOOST_CONCEPT_ASSERT((ReadWritePropertyMapConcept<Bags, Decomposition_vertex>));
        //Graph_vertex and Bag_inside_type must be a same type
        BOOST_STATIC_ASSERT((is_same<Graph_vertex, Bag_inside_type>::value));

        subgraph<Graph> sg;
        copy_graph(g, sg);
        std::vector<subgraph<Graph>*> components = detail::get_components(sg);

        bool first_component = true;
        Decomposition_vertex prev_root;
        for (subgraph<Graph>* curr : components) {
            //create w set of vertices
            std::set<Graph_vertex> w;
            for (const auto v : make_iterator_range(vertices(*curr))) w.insert((*curr).local_to_global(v));
            tuple<Decomposition_vertex, bool> res = detail::decompose(*curr, d, bags, w, std::set<Graph_vertex>(), k);
            if (!get<1>(res)) return false;
            if (first_component) first_component = false;
            //if it's not a first component, we connect current root with root of previous component
            else add_edge(get<0>(res), prev_root, d);
            prev_root = get<0>(res);
        }

        return true;
    }
} // namespace boost
#endif // BOOST_GRAPH_TREE_DECOMPOSITION_HPP
