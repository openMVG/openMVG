
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "testing/testing.h"

#include "dynamic_bitset.hpp"

TEST(DYNAMIC_BITSET, InitAndReset_64)
{
  using namespace stl;

  const int nbBits = 64;
  dynamic_bitset mybitset(nbBits);

  // Check that there is nbBits bit stored
  EXPECT_EQ(64, mybitset.size());
  // Check that there is just the necessary count of BlockType allocated for storage
  EXPECT_EQ(64/dynamic_bitset::bits_per_block, mybitset.num_blocks());
  
  // Set some bits to 1
  for (int i = 0; i < mybitset.size(); i+=2)
    mybitset[i] = true;

  // Check that some bits have been correctly set to 1
  for (int i = 0; i < mybitset.size(); ++i)
  {
    EXPECT_EQ(!(i%2), mybitset[i]);
  }
  
  // Reset the value to 0
  mybitset.reset();
  for (int i = 0; i < mybitset.size(); ++i)
  {
    EXPECT_EQ(false, mybitset[i]);
  }
}

// Create a dynamic_bitset that is shorter than the internal used bit container
TEST(DYNAMIC_BITSET, InitAndReset_4)
{
  using namespace stl;

  const int nbBits = 4;
  dynamic_bitset mybitset(nbBits);
  
  EXPECT_EQ(4, mybitset.size());
  EXPECT_EQ(1, mybitset.num_blocks());

  // Set some bits to 1
  for (int i = 0; i < mybitset.size(); i+=2)
    mybitset[i] = true;

  // Check that some bits have been correctly set to 1
  for (int i = 0; i < mybitset.size(); ++i)
  {
    EXPECT_EQ(!(i%2), mybitset[i]);
  }
  
  // Reset the value to 0
  mybitset.reset();
  for (int i = 0; i < mybitset.size(); ++i)
  {
    EXPECT_EQ(false, mybitset[i]);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

