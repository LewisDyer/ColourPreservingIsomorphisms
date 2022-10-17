#ifndef BOOST_GRAPH_NICE_TREE_DECOMPOSITION_HPP
#define BOOST_GRAPH_NICE_TREE_DECOMPOSITION_HPP

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

/**
 * @file nice_tree_decomposition.hpp
 *
 * The function \p nice_tree_decomposition() constructs a nice tree decomposition based on decomposition \p d of same width and stores it in provided graph \p nice_d. Bags for each node will be stored in property map \p nice_bags.
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
         * Retrieves neighbour vertices of target_v that haven't been processed yet
         * @tparam Decomposition    Type of the decomposition
         * @param d                 Decomposition in which neighbours are looked up
         * @param target_v          Vertex which neighbours will be retrieved
         * @param visitor_map       Map of vertices, where bool=true means that vertex has been already processed
         * @return List of neighbours of target_v in dš
         */
        template <class Decomposition>
        std::list<typename graph_traits<Decomposition>::vertex_descriptor>
        get_neighbour_vertices(Decomposition & d,
                               typename graph_traits<Decomposition>::vertex_descriptor target_v,
                               std::map<typename graph_traits<Decomposition>::vertex_descriptor, bool> & visitor_map) {
            //vertex descriptor that is used to address vertices in decomposition graph
            using Decomposition_vertex = typename graph_traits<Decomposition>::vertex_descriptor;

            std::list<Decomposition_vertex> neighbours;
            auto adj = adjacent_vertices(target_v, d);
            for (const Decomposition_vertex neighbour : make_iterator_range(adj)) {
                if (!visitor_map[neighbour]) neighbours.push_back(neighbour);
            }

            return neighbours;
        }

        /**
         * Recursive function that constructs nice tree decomposition from general decomposition
         * @tparam Decomposition        Type of the decomposition
         * @tparam Bags                 Type of the bags of Decomposition
         * @param d                     General decomposition provided by the user
         * @param visitor_map           Map of vertices, where bool=true means that vertex has been already processed
         * @param nice_d                Nice tree decomposition that's being constructed
         * @param bags                  Bags of d
         * @param nice_bags             Bags of nice_d
         * @param d_head_vertex         Vertex of d that's being currently processed
         * @param nice_d_tail_vertex    Vertex of nice_d that's being currently processed
         */
        template <class Decomposition, class Bags>
        void nice_decomposition_rec(Decomposition & d,
                                    std::map<typename graph_traits<Decomposition>::vertex_descriptor, bool> & visitor_map,
                                    Decomposition & nice_d,
                                    Bags & bags,
                                    Bags & nice_bags,
                                    typename graph_traits<Decomposition>::vertex_descriptor d_head_vertex,
                                    typename graph_traits<Decomposition>::vertex_descriptor nice_d_tail_vertex) {
            //vertex descriptor that is used to address vertices in decomposition graph
            using Decomposition_vertex = typename graph_traits<Decomposition>::vertex_descriptor;
            using Vertex_descriptor_list = std::list<Decomposition_vertex>;

            //marking current d_head_vertex as visited
            visitor_map[d_head_vertex] = true;

            //checking if there is vertex in (nice_d_tail_vertex \ d_head_vertex)
            for (const Decomposition_vertex curr_v : nice_bags[nice_d_tail_vertex]) {
                if (std::find(bags[d_head_vertex].begin(), bags[d_head_vertex].end(), curr_v) == bags[d_head_vertex].end()) {
                    //creating new nice decomposition node
                    Decomposition_vertex new_v = add_vertex(nice_d);
                    add_edge(nice_d_tail_vertex, new_v, nice_d);
                    //creating bag without erasing curr_v
                    for (const Decomposition_vertex copy_v : nice_bags[nice_d_tail_vertex]) {
                        if (copy_v != curr_v) nice_bags[new_v].insert(nice_bags[new_v].begin(), copy_v);
                    }
                    nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, d_head_vertex, new_v);
                    return;
                }
            }

            //checking if there is vertex in (d_head_vertex \ nice_d_tail_vertex)
            for (const Decomposition_vertex curr_v : bags[d_head_vertex]) {
                if (std::find(nice_bags[nice_d_tail_vertex].begin(), nice_bags[nice_d_tail_vertex].end(), curr_v) == nice_bags[nice_d_tail_vertex].end()) {
                    //creating new nice decomposition node
                    Decomposition_vertex new_v = add_vertex(nice_d);
                    add_edge(nice_d_tail_vertex, new_v, nice_d);
                    //creating bag and inserting curr_v
                    for (const Decomposition_vertex copy_v : nice_bags[nice_d_tail_vertex]) nice_bags[new_v].insert(nice_bags[new_v].begin(), copy_v);
                    nice_bags[new_v].insert(nice_bags[new_v].begin(), curr_v);
                    nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, d_head_vertex, new_v);
                    return;
                }
            }

            //if call reached here, it means nice_d_tail_vertex and d_head_vertex are equal
            //finding neighbours of d_head_vertex
            Vertex_descriptor_list neighbours = get_neighbour_vertices(d, d_head_vertex, visitor_map);
            if (neighbours.size() == 0) {
                Vertex_descriptor_list remaining_v;
                for (const Decomposition_vertex v : bags[d_head_vertex]) remaining_v.push_back(v);
                //inserting new nodes to nice decomposition until remaining_v is empty
                while (remaining_v.size()) {
                    remaining_v.pop_front();
                    //creating new nice decomposition node
                    Decomposition_vertex new_v = add_vertex(nice_d);
                    add_edge(nice_d_tail_vertex, new_v, nice_d);
                    //adding elements of remaining_v to newly created bag of nice decomposition
                    for (const Decomposition_vertex v :remaining_v) nice_bags[new_v].insert(nice_bags[new_v].begin(), v);
                    nice_d_tail_vertex = new_v;
                }
            }
            else if (neighbours.size() == 1) nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, neighbours.front(), nice_d_tail_vertex);
            else if (neighbours.size() == 2) {
                for (const Decomposition_vertex neighbour : neighbours) {
                    //creating new nice decomposition node
                    Decomposition_vertex new_v = add_vertex(nice_d);
                    add_edge(nice_d_tail_vertex, new_v, nice_d);
                    //adding remaining_v to newly created bag of nice decomposition
                    for (const Decomposition_vertex copy_v : bags[d_head_vertex]) nice_bags[new_v].insert(nice_bags[new_v].begin(), copy_v);
                    nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, neighbour, new_v);
                }
            }
            else {
                Decomposition_vertex first_neighbour = neighbours.front();
                neighbours.pop_front();
                //creating new nice decomposition node
                Decomposition_vertex new_v = add_vertex(nice_d);
                add_edge(nice_d_tail_vertex, new_v, nice_d);
                //adding remaining_v to newly created bag of nice decomposition
                for (const Decomposition_vertex copy_v : bags[d_head_vertex]) nice_bags[new_v].insert(nice_bags[new_v].begin(), copy_v);
                //calling recursion for first neighbour
                nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, first_neighbour, new_v);

                //creating parent node for other neighbours
                new_v = add_vertex(nice_d);
                add_edge(nice_d_tail_vertex, new_v, nice_d);
                //adding remaining_v to newly created bag of nice decomposition
                for (const Decomposition_vertex copy_v : bags[d_head_vertex]) nice_bags[new_v].insert(nice_bags[new_v].begin(), copy_v);
                nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, d_head_vertex, new_v);
            }
        }

        /**
         * Inits the recursion that constructs the nice tree decomposition
         * @tparam Decomposition    Type of the decomposition
         * @tparam Bags             Type of the bags of Decomposition
         * @param d                 General decomposition provided by the user
         * @param nice_d            Nice tree decomposition that's being constructed
         * @param bags              Bags of d
         * @param nice_bags         Bags of nice_d
         * @return Root vertex of nice_d
         */
        template <class Decomposition, class Bags>
        typename graph_traits<Decomposition>::vertex_descriptor
        nice_tree_decomposition_dispatch(Decomposition & d,
                                         Decomposition & nice_d,
                                         Bags & bags,
                                         Bags & nice_bags) {
            //vertex descriptor that is used to address vertices in decomposition graph
            using Decomposition_vertex = typename graph_traits<Decomposition>::vertex_descriptor;
            //visitor map of decomposition d
            using Visitor_map = std::map<Decomposition_vertex, bool>;

            //creating visitor map
            Visitor_map visitor_map;
            for (const Decomposition_vertex v : make_iterator_range(vertices(d))) visitor_map[v] = false;
            //creating root of nice decomposition
            Decomposition_vertex root = add_vertex(nice_d);
            //getting first vertex of d
            Decomposition_vertex first_v = *(vertices(d).first);
            //calling recursive function that builds nice decomposition
            nice_decomposition_rec(d, visitor_map, nice_d, bags, nice_bags, first_v, root);

            return root;
        }

        /**
         * Checks decomposition provided by the user for cycles - if cycle is found, std::invalid_argument is thrown
         * @tparam Decomposition    Type of the decomposition
         * @param d                 General decomposition provided by the user
         * @throws std::invalid_argument If d contains a cycle
         */
        template <class Decomposition>
        void cycle_check(Decomposition & d) {
            //vertex descriptor that is used to address edges in decomposition graph
            using Decomposition_edge_descriptor = typename graph_traits<Decomposition>::edge_descriptor;

            //visitor class that detects cycle and throws exception
            class dfs_cycle_visitor : public dfs_visitor<> {
                public:
                    void tree_edge(Decomposition_edge_descriptor e, const Decomposition & /* g */) {
                        m_tree_edges.insert(e);
                    }

                    void back_edge(Decomposition_edge_descriptor e, const Decomposition & /* g */) const {
                        if (m_tree_edges.find(e) == m_tree_edges.end()) throw std::invalid_argument("Invalid nice tree decomposition - decomposition contains a cycle");
                    }
                private:
                    std::set<Decomposition_edge_descriptor> m_tree_edges;
            };

            if (num_vertices(d) == 0) return;
            dfs_cycle_visitor vis{};
            depth_first_search(d, visitor(vis));
        }
    } // namespace detail
    #endif // DOXYGEN_SHOULD_SKIP_THIS

    /**
     * @tparam      Decomposition   Type of the decomposition.
     * @tparam      Bags            Type of the property map representing bags of \p Decomposition.
     * @param[in]   d               An undirected graph containing decomposition for which will be the nice decomposition constructed. The graph type must be a model of VertexListGraph.
     * @param[in]   bags            The bags property map representing bags of \p d. The type must be a model of a mutable Readable Property Map. The key type of the map must be a vertex descriptor of \p d. The value type of the map must be a model of Mutable_Container. It also must implement function insert([position/hint iterator], [element]).
     * @param[out]  nice_d          An undirected graph in which the nice decomposition will be stored at the end of algorithm.
     * @param[out]  nice_bags       The bags property map representing bags of \p nice_d. The type must be a model of a mutable Readable Property Map. The key type of the map must be a vertex descriptor of \p d. The value type of the map must be a model of Mutable_Container. It also must implement function insert([position/hint iterator], [element]).
     * @throws std::invalid_argument If \p d contains a cycle
     * @return Vertex descriptor of root of the nice decomposition stored in \p nice_d.
     */
    template <class Decomposition, class Bags>
    typename graph_traits<Decomposition>::vertex_descriptor nice_tree_decomposition(Decomposition & d, Bags & bags, Decomposition & nice_d, Bags & nice_bags) {
        //vertex descriptor that is used to address vertices in decomposition graph
        using Decomposition_vertex = typename graph_traits<Decomposition>::vertex_descriptor;
        //type of value (container) of PropertyMap
        using Bag_value_type = typename property_traits<Bags>::value_type;

        //Decomposition must be model of VertexListGraphConcept
        BOOST_CONCEPT_ASSERT((VertexListGraphConcept<Decomposition>));
        //each bag of Bags container must be a model of Mutable_Container
        BOOST_CONCEPT_ASSERT((detail::InsertCollection<Bag_value_type>));
        //Bags must be Read Write Property Map (with Decomposition_vertex as a key)
        BOOST_CONCEPT_ASSERT((ReadWritePropertyMapConcept<Bags, Decomposition_vertex>));

        detail::cycle_check(d);

        return detail::nice_tree_decomposition_dispatch(d, nice_d, bags, nice_bags);
    }
} // namespace boost
#endif // BOOST_GRAPH_NICE_TREE_DECOMPOSITION_HPP
