// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/triplet_finder.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <iostream>
#include <vector>

using namespace openMVG::graph;

using Pairs = std::vector<std::pair<int,int>>;

TEST(TripletFinder, test_no_triplet) {

  // a - b - c
  const int a = 0, b = 1, c = 2;
  const Pairs pairs = {{a, b}, {b, c}};

  std::vector<Triplet> vec_triplets;
  EXPECT_FALSE(ListTriplets(pairs, vec_triplets));
  EXPECT_TRUE(vec_triplets.empty());
}

TEST(TripletFinder, test_one_triplet) {

  {

    //
    // a_b
    // |/
    // c
    const int a = 0, b = 1, c = 2;
    const Pairs pairs = {{a, b}, {a, c}, {b, c}};

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(1, vec_triplets.size());
    //Check the cycle values
    EXPECT_EQ(0, vec_triplets[0].i);
    EXPECT_EQ(1, vec_triplets[0].j);
    EXPECT_EQ(2, vec_triplets[0].k);
  }

  {
    //
    // a_b__c
    //    |/
    //    d
    const int a = 0, b = 1, c = 2, d = 3;
    const Pairs pairs = {{a, b}, {b, c}, {b, d}, {c, d}};

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(1, vec_triplets.size());
    //Check the cycle values
    EXPECT_EQ(1,vec_triplets[0].i);
    EXPECT_EQ(2,vec_triplets[0].j);
    EXPECT_EQ(3,vec_triplets[0].k);
  }
}

TEST(TripletFinder, test_two_triplet) {

  {
    //
    // a__b
    // |\ |
    // | \|
    // c--d
    const int a = 0, b = 1, c = 2, d = 3;
    const Pairs pairs = {{a, b}, {a, c}, {a, d}, {c, d}, {b,d}};

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(2, vec_triplets.size());
  }

  {
    //
    // a   c
    // |\ /|
    // | b |
    // |/ \|
    // d   e
    const int a = 0, b = 1, c = 2, d = 3, e = 4;
    const Pairs pairs = {{a, b}, {b,c}, {c,e}, {e,b}, {b,d}, {d,a}};

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(2, vec_triplets.size());
  }

    {
    //
    // a      c
    // |\    /|
    // | b--f |
    // |/    \|
    // d      e
    const int a = 0, b = 1, c = 2, d = 3, e = 4, f = 5;
    const Pairs pairs = {{a,b}, {b,f}, {f,c}, {c,e}, {e,f}, {f,b}, {b,d}, {d,a}};

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(2, vec_triplets.size());
  }
}

TEST(TripletFinder, test_three_triplet) {

  {
    //
    // a   b
    // |\ /|
    // c-d-e
    // |/
    // f
    const int a = 0, b = 1, c = 2, d = 3, e = 4, f = 5;
    const Pairs pairs = { {a,c}, {a,d}, {c,d}, {c,f}, {f,d}, {d,b}, {b,e}, {e,d} };

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(3, vec_triplets.size());
  }

  {
    //
    // a        b--g--h
    // | \    / |   \/
    // |  d--e  |    i
    // | /    \ |
    // c        f

    const int a = 0, b = 1, c = 2, d = 3, e = 4, f = 5, g = 6, h = 7, i = 8;
    const Pairs pairs = {
      {a,c}, {a,d}, {d,c},
      {d,e},
      {e,b}, {e,f}, {b,f},
      {b,g},
      {g,h}, {h,i}, {i,g} };

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(3, vec_triplets.size());
  }

  {
    //
    // a---b
    // |\  |\
    // | \ | \
    // |  \|  \
    // c---d---e
    //
    const int a = 0, b = 1, c = 2, d = 3, e = 4;
    const Pairs pairs = { {a,b}, {b,d}, {d,c}, {c,a}, {a,d}, {b,e}, {d,e} };

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(3, vec_triplets.size());
  }

}

TEST(TripletFinder, test_for_triplet) {

 {
    //
    // a__b
    // |\/|
    // |/\|
    // c--d
    const int a = 0, b = 1, c = 2, d = 3;
    const Pairs pairs = { {a,b}, {a,c}, {a,d}, {c,d}, {b,d}, {c,b} };

    std::vector<Triplet> vec_triplets;
    EXPECT_TRUE(ListTriplets(pairs, vec_triplets));
    EXPECT_EQ(4, vec_triplets.size());
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
