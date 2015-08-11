// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "testing/testing.h"
#include "openMVG/matching/metric.hpp"
#include <iostream>
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
  
  Hamming< unsigned char > metricHammingUchar;
  
  EXPECT_EQ(8, metricHammingUchar(reinterpret_cast<unsigned char *>(&a),reinterpret_cast<unsigned char *>(&b),1));
  EXPECT_EQ(0, metricHammingUchar(reinterpret_cast<unsigned char *>(&a),reinterpret_cast<unsigned char *>(&a),1));
  EXPECT_EQ(2, metricHammingUchar(reinterpret_cast<unsigned char *>(&a),reinterpret_cast<unsigned char *>(&c),1));
}

TEST(Metric, HAMMING_BITSET_RAW_MEMORY_64BITS)
{
  const int COUNT = 4;
  std::bitset<64> tab[COUNT];
  // Zeros
  for(int i = 0; i < 64; ++i) {  tab[0][i] = 0;  }
  // 0101 ...
  for(int i = 0; i < 64; ++i) {  tab[1][i] = i%2 == 0;  }
  // 00110011...
  for(int i = 0; i < 64; ++i) {  tab[2][i] = (i/2)%2 == 0;  }
  // 000111000111...
  for(int i = 0; i < 64; ++i) {  tab[3][i] = (i/3)%2 == 0;  }

  // ground truth hamming distances between bit array
  const double gtDist[] =
  {0, 32, 32, 33, 32,
   0, 32, 21, 32, 32,
   0, 31, 33, 21, 31, 0};

  HammingBitSet<std::bitset<8> > metricHammingBitSet;
  Hamming< unsigned char > metricHamming;
  size_t cpt = 0;
  for (size_t i = 0; i < COUNT; ++i)
  {
    for (size_t j = 0; j < COUNT; ++j)
    {
      EXPECT_EQ(gtDist[cpt], metricHammingBitSet(&tab[i],&tab[j], 1));
      EXPECT_EQ(gtDist[cpt], metricHamming((uint64_t*)&tab[i],(uint64_t*)&tab[j], sizeof(uint64_t)));
      //Check distance symmetry
      EXPECT_EQ(gtDist[cpt], metricHammingBitSet(&tab[j],&tab[i], 1));
      EXPECT_EQ(gtDist[cpt], metricHamming((uint64_t*)&tab[j],(uint64_t*)&tab[i], sizeof(uint64_t)));
      ++cpt;
    }
  }
}

TEST(Metric, HAMMING_BITSET_RAW_MEMORY_32BITS)
{
  const int COUNT = 4;
  std::bitset<32> tab[COUNT];
  // Zeros
  for(int i = 0; i < 32; ++i) {  tab[0][i] = 0;  }
  // 0101 ...
  for(int i = 0; i < 32; ++i) {  tab[1][i] = i%2 == 0;  }
  // 00110011...
  for(int i = 0; i < 32; ++i) {  tab[2][i] = (i/2)%2 == 0;  }
  // 000111000111...
  for(int i = 0; i < 32; ++i) {  tab[3][i] = (i/3)%2 == 0;  }

  // ground truth hamming distances between bit array
  const double gtDist[] =
  {0, 16, 16, 17, 16,
  0, 16, 11, 16, 16,
  0, 17, 17, 11, 17, 0};

  HammingBitSet<std::bitset<8> > metricHammingBitSet;
  Hamming< unsigned char > metricHamming;
  size_t cpt = 0;
  for (size_t i = 0; i < COUNT; ++i)
  {
    for (size_t j = 0; j < COUNT; ++j)
    {
      EXPECT_EQ(gtDist[cpt], metricHammingBitSet(&tab[i],&tab[j], 1));
      EXPECT_EQ(gtDist[cpt], metricHamming((uint32_t*)&tab[i],(uint32_t*)&tab[j], sizeof(uint32_t)));
      //Check distance symmetry
      EXPECT_EQ(gtDist[cpt], metricHammingBitSet(&tab[j],&tab[i], 1));
      EXPECT_EQ(gtDist[cpt], metricHamming((uint32_t*)&tab[j],(uint32_t*)&tab[i], sizeof(uint32_t)));
      ++cpt;
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
