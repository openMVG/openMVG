// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "testing/testing.h"
#include "openMVG/matching/metric.hpp"
#include "third_party/flann/src/cpp/flann/algorithms/dist.h"
#include <iostream>
#include <bitset>
#include <string>
using namespace std;

using namespace openMVG;
using namespace matching;

template<typename Metric>
typename Metric::ResultType DistanceT()
{
  typename Metric::ElementType array1[] = {0, 1, 2, 3, 4, 5, 6, 7};
  typename Metric::ElementType array2[] = {7, 6, 5, 4, 3, 2, 1, 0};
  Metric metric;
  return metric(array1, array2, 8);
}

TEST(Metric, L2_Simple)
{
  EXPECT_EQ(168, DistanceT<L2_Simple<unsigned char> >());
  EXPECT_EQ(168, DistanceT<L2_Simple<short> >());
  EXPECT_EQ(168, DistanceT<L2_Simple<int> >());
  EXPECT_EQ(168, DistanceT<L2_Simple<float> >());
  EXPECT_EQ(168, DistanceT<L2_Simple<double> >());
}

TEST(Metric, L2_Vectorized)
{
  EXPECT_EQ(168, DistanceT<L2_Vectorized<unsigned char> >());
  EXPECT_EQ(168, DistanceT<L2_Vectorized<short> >());
  EXPECT_EQ(168, DistanceT<L2_Vectorized<int> >());
  EXPECT_EQ(168, DistanceT<L2_Vectorized<float> >());
  EXPECT_EQ(168, DistanceT<L2_Vectorized<double> >());
}

TEST(Metric, HAMMING_BITSET)
{
  std::bitset<8> a(std::string("01010101"));
  std::bitset<8> b(std::string("10101010"));
  std::bitset<8> c(std::string("11010100"));

  HammingBitSet<std::bitset<8> > metricHamming;
  EXPECT_EQ(8, metricHamming(&a,&b,1));
  EXPECT_EQ(0, metricHamming(&a,&a,1));
  EXPECT_EQ(2, metricHamming(&a,&c,1));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
