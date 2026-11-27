// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/matching/metric.hpp"
#if defined(OPENMVG_USE_AVX2) || defined(OPENMVG_USE_AVX)
#include "openMVG/system/cpu_instruction_set.hpp"
#endif

#include "testing/testing.h"

#include <array>
#include <iostream>

using namespace openMVG;
using namespace matching;

template<typename Metric>
typename Metric::ResultType DistanceT()
{
  typename Metric::ElementType array1[] = {0, 1, 2, 3, 4, 5, 6, 7};
  typename Metric::ElementType array2[] = {7, 6, 5, 4, 3, 2, 1, 0};
  const Metric metric{};
  return metric(array1, array2, 8);
}

TEST(Metric, HAMMING_UNALIGNED_RAW_MEMORY)
{
  alignas(8) std::array<uint8_t, 64 + 7> a, b;
  for (size_t i = 0; i < a.size(); ++i)
  {
    a[i] = static_cast<uint8_t>(37 * i + 19);
    b[i] = static_cast<uint8_t>(17 * i + 203);
  }

  const Hamming<uint8_t> metric{};
  // Exercise the 64-bit, 32-bit and byte paths at every word alignment.
  for (size_t size : {0, 1, 3, 4, 7, 8, 12, 28, 32, 64})
  {
    for (size_t offset_a = 0; offset_a < 8; ++offset_a)
    {
      for (size_t offset_b = 0; offset_b < 8; ++offset_b)
      {
        unsigned int expected = 0;
        for (size_t i = 0; i < size; ++i)
          expected += std::bitset<8>(a[offset_a + i] ^ b[offset_b + i]).count();

        EXPECT_EQ(expected, metric(a.data() + offset_a, b.data() + offset_b, size));
        EXPECT_EQ(expected, metric(b.data() + offset_b, a.data() + offset_a, size));
        EXPECT_EQ(0, metric(a.data() + offset_a, a.data() + offset_a, size));
      }
    }
  }
}

TEST(Metric, L2DIM128_UNALIGNED_RAW_MEMORY)
{
  alignas(32) std::array<uint8_t, 128 + 31> a, b;
  for (size_t i = 0; i < a.size(); ++i)
  {
    a[i] = static_cast<uint8_t>(37 * i + 19);
    b[i] = static_cast<uint8_t>(17 * i + 203);
  }

  const L2<uint8_t> metric{};
  for (size_t offset_a = 0; offset_a < 32; ++offset_a)
  {
    for (size_t offset_b = 0; offset_b < 32; ++offset_b)
    {
      int expected = 0;
      for (size_t i = 0; i < 128; ++i)
      {
        const int difference = int(a[offset_a + i]) - int(b[offset_b + i]);
        expected += difference * difference;
      }
      EXPECT_EQ(expected, metric(a.data() + offset_a, b.data() + offset_b, 128));
      EXPECT_EQ(expected, metric(b.data() + offset_b, a.data() + offset_a, 128));
      EXPECT_EQ(0, metric(a.data() + offset_a, a.data() + offset_a, 128));
    }
  }
}

TEST(Metric, L2)
{
  EXPECT_EQ(168, DistanceT<L2<uint8_t>>());
  EXPECT_EQ(168, DistanceT<L2<uint16_t>>());
  EXPECT_EQ(168, DistanceT<L2<int16_t>>());
  EXPECT_EQ(168, DistanceT<L2<int32_t>>());
  EXPECT_EQ(168.f, DistanceT<L2<float>>());
  EXPECT_EQ(168.0, DistanceT<L2<double>>());
}

TEST(Metric, HAMMING_BITSET)
{
  const std::bitset<8>
    a("01010101"),
    b("10101010"),
    c("11010100");

  const HammingBitSet<std::bitset<8>> metricHamming{};
  EXPECT_EQ(8, metricHamming(&a,&b,1));
  EXPECT_EQ(0, metricHamming(&a,&a,1));
  EXPECT_EQ(2, metricHamming(&a,&c,1));

  const Hamming<unsigned char> metricHammingUchar{};

  EXPECT_EQ(8, metricHammingUchar(reinterpret_cast<const uint8_t *>(&a),reinterpret_cast<const uint8_t *>(&b),1));
  EXPECT_EQ(0, metricHammingUchar(reinterpret_cast<const uint8_t *>(&a),reinterpret_cast<const uint8_t *>(&a),1));
  EXPECT_EQ(2, metricHammingUchar(reinterpret_cast<const uint8_t *>(&a),reinterpret_cast<const uint8_t *>(&c),1));
}

TEST(Metric, HAMMING_BITSET_RAW_MEMORY_64BITS)
{
  const int COUNT = 4;
  std::bitset<64> tab[COUNT];
  // Zeros
  for (int i = 0; i < 64; ++i) {  tab[0][i] = 0;  }
  // 0101 ...
  for (int i = 0; i < 64; ++i) {  tab[1][i] = i%2 == 0;  }
  // 00110011...
  for (int i = 0; i < 64; ++i) {  tab[2][i] = (i/2)%2 == 0;  }
  // 000111000111...
  for (int i = 0; i < 64; ++i) {  tab[3][i] = (i/3)%2 == 0;  }

  // ground truth hamming distances between bit array
  const double gtDist[] =
    {0, 32, 32, 33, 32,
     0, 32, 21, 32, 32,
     0, 31, 33, 21, 31, 0};

  HammingBitSet<std::bitset<8>> metricHammingBitSet;
  const Hamming<unsigned char> metricHamming{};
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
  for (int i = 0; i < 32; ++i) {  tab[0][i] = 0;  }
  // 0101 ...
  for (int i = 0; i < 32; ++i) {  tab[1][i] = i%2 == 0;  }
  // 00110011...
  for (int i = 0; i < 32; ++i) {  tab[2][i] = (i/2)%2 == 0;  }
  // 000111000111...
  for (int i = 0; i < 32; ++i) {  tab[3][i] = (i/3)%2 == 0;  }

  // ground truth hamming distances between bit array
  const double gtDist[] =
    {0, 16, 16, 17, 16,
    0, 16, 11, 16, 16,
    0, 17, 17, 11, 17, 0};

  HammingBitSet<std::bitset<8>> metricHammingBitSet;
  const Hamming<unsigned char> metricHamming{};
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

TEST(Metric, L2DIM128)
{
  // Test SIFT like descriptor (uint8_t)
  {
    using VecUC128 = Eigen::Matrix<uint8_t, 128, 1>;
    const VecUC128 a = VecUC128::Random();
    const VecUC128 b = VecUC128::Random();
    const unsigned int GTL2 = (a.cast<int>()-b.cast<int>()).squaredNorm();
    const L2<uint8_t> metricL2{};
    EXPECT_EQ(GTL2, metricL2(a.data(), b.data(), 128));
    #ifdef OPENMVG_USE_AVX2
      openMVG::system::CpuInstructionSet cpu_instruction_set;
      EXPECT_TRUE(cpu_instruction_set.supportAVX2());
      EXPECT_EQ(GTL2, L2_AVX2(a.data(), b.data(), 128));
    #endif
  }

  // Test SIFT like descriptor (float)
  {
    using VecF128 = Eigen::Matrix<float, 128, 1>;
    const VecF128 a = VecF128::Random();
    const VecF128 b = VecF128::Random();
    const double GTL2 = (a-b).squaredNorm();
    const L2<float> metricL2{};
    EXPECT_NEAR(GTL2, metricL2(a.data(), b.data(), 128), 1e-4);
    #ifdef OPENMVG_USE_AVX
      openMVG::system::CpuInstructionSet cpu_instruction_set;
      EXPECT_TRUE(cpu_instruction_set.supportAVX());
      EXPECT_NEAR(GTL2, L2_AVX(a.data(), b.data(), 128), 1e-4);
    #endif
  }
}

TEST(Metric, L1DIM128) {
    using VecUC128 = Eigen::Matrix<uint8_t, 128, 1>;
    const VecUC128 a = VecUC128::Random();
    const VecUC128 b = VecUC128::Random();
    const unsigned int GTL1 = (a.cast<int>()-b.cast<int>()).lpNorm<1>();
    const L1<uint8_t> metricL1{};
    EXPECT_NEAR(GTL1, metricL1(a.data(), b.data(), 128), 1e-4);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
