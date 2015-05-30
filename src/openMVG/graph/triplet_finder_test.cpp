// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/triplet_finder.hpp"
using namespace openMVG::graph;

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <iostream>
#include <vector>

TEST(TripletFinder, test_no_triplet) {

  typedef lemon::ListGraph Graph;

  // a_b_c
  Graph ga;
  Graph::Node a = ga.addNode(), b = ga.addNode(), c = ga.addNode();
  ga.addEdge(a,b);
  ga.addEdge(b,c);

  std::vector< Triplet > vec_triplets;

  EXPECT_FALSE(List_Triplets(ga, vec_triplets));
  EXPECT_TRUE(vec_triplets.empty());
}

TEST(TripletFinder, test_one_triplet) {

  typedef lemon::ListGraph Graph;

  {
    //
    // a_b
    // |/
    // c
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(), c = ga.addNode();
    ga.addEdge(a,b);
    ga.addEdge(a,c);
    ga.addEdge(b,c);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(1, vec_triplets.size());
    //Check the cycle values
    EXPECT_EQ(0,vec_triplets[0].i);
    EXPECT_EQ(1,vec_triplets[0].j);
    EXPECT_EQ(2,vec_triplets[0].k);
  }

  {
    //
    // a_b__c
    //    |/
    //    d
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode();
    ga.addEdge(a,b);
    ga.addEdge(b,c);
    ga.addEdge(b,d);
    ga.addEdge(c,d);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(1, vec_triplets.size());
    //Check the cycle values
    EXPECT_EQ(1,vec_triplets[0].i);
    EXPECT_EQ(2,vec_triplets[0].j);
    EXPECT_EQ(3,vec_triplets[0].k);
  }
}

TEST(TripletFinder, test_two_triplet) {

  typedef lemon::ListGraph Graph;

  {
    //
    // a__b
    // |\ |
    // | \|
    // c--d
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode();

    ga.addEdge(a,b);
    ga.addEdge(a,c);
    ga.addEdge(a,d);
    ga.addEdge(c,d);
    ga.addEdge(b,d);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(2, vec_triplets.size());
  }

  {
    //
    // a   c
    // |\ /|
    // | b |
    // |/ \|
    // d   e
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode(),
      e = ga.addNode();

    ga.addEdge(a,b);
    ga.addEdge(b,c);
    ga.addEdge(c,e);
    ga.addEdge(e,b);
    ga.addEdge(b,d);
    ga.addEdge(d,a);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(2, vec_triplets.size());
  }

    {
    //
    // a      c
    // |\    /|
    // | b--f |
    // |/    \|
    // d      e
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode(),
      e = ga.addNode(), f = ga.addNode();

    ga.addEdge(a,b);
    ga.addEdge(b,f);
    ga.addEdge(f,c);
    ga.addEdge(c,e);
    ga.addEdge(e,f);
    ga.addEdge(f,b);
    ga.addEdge(b,d);
    ga.addEdge(d,a);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(2, vec_triplets.size());
  }
}


TEST(TripletFinder, test_three_triplet) {

  typedef lemon::ListGraph Graph;

  {
    //
    // a   b
    // |\ /|
    // c-d-e
    // |/
    // f
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode(),
      e = ga.addNode(), f = ga.addNode();

    ga.addEdge(a,c);
    ga.addEdge(a,d);
    ga.addEdge(c,d);
    ga.addEdge(c,f);
    ga.addEdge(f,d);
    ga.addEdge(d,b);
    ga.addEdge(b,e);
    ga.addEdge(e,d);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(3, vec_triplets.size());
  }

  {
    //
    // a        b--g--h
    // | \    / |   \/
    // |  d--e  |    i
    // | /    \ |
    // c        f
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode(),
      e = ga.addNode(), f= ga.addNode(),
      g = ga.addNode(), h = ga.addNode(),
      i = ga.addNode();

    ga.addEdge(a,c);
    ga.addEdge(a,d);
    ga.addEdge(d,c);
    ga.addEdge(d,e);
    ga.addEdge(e,b);
    ga.addEdge(e,f);
    ga.addEdge(b,f);
    ga.addEdge(b,g);
    ga.addEdge(g,h);
    ga.addEdge(h,i);
    ga.addEdge(i,g);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(3, vec_triplets.size());
  }

    {
    //
    // a---b
    // |\  |\
    // | \ | \
    // |  \|  \
    // c---d---e
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode(),
      e = ga.addNode();

    ga.addEdge(a,b);
    ga.addEdge(b,d);
    ga.addEdge(d,c);
    ga.addEdge(c,a);
    ga.addEdge(a,d);
    ga.addEdge(b,e);
    ga.addEdge(d,e);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(3, vec_triplets.size());
  }
}

TEST(TripletFinder, test_for_triplet) {

  typedef lemon::ListGraph Graph;
 {
    //
    // a__b
    // |\/|
    // |/\|
    // c--d
    Graph ga;
    Graph::Node a = ga.addNode(), b = ga.addNode(),
      c = ga.addNode(), d = ga.addNode();

    ga.addEdge(a,b);
    ga.addEdge(a,c);
    ga.addEdge(a,d);
    ga.addEdge(c,d);
    ga.addEdge(b,d);
    ga.addEdge(c,b);

    std::vector< Triplet > vec_triplets;
    EXPECT_TRUE(List_Triplets(ga, vec_triplets));
    EXPECT_EQ(4, vec_triplets.size());
  }
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
