#include <gtest/gtest.h>
#include "../count_colourful_isomorphisms.h"
// using namespace boost;
// using Graph = adjacency_list<listS, vecS, undirectedS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;
// using DiGraph = adjacency_list<listS, vecS, bidirectionalS, property<vertex_name_t, int, property<vertex_index_t, size_t>>, property<edge_index_t, int>>;



class TreeCPITest : public testing::Test {
  protected:
    void SetUp() override {
      // some testing here
      s = star<Graph>(4);
      g_star = makeRooted(s, 0);
      for(int i=0; i < boost::num_vertices(g_star); i++) {col_star[i] = i;}

      p = path<Graph>(10);
      for(int i=0; i < boost::num_vertices(p); i++) {
        starInPath[i] = i % boost::num_vertices(g_star);
      }
      

    }

    Graph s, p;
    DiGraph g_star;
    Colours col_star, starInPath;
    
};

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

TEST_F(TreeCPITest, StarNotFoundInPath) {
  EXPECT_EQ(tree_count(g_star, 0, col_star, p, starInPath), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}