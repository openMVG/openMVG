// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/features.hpp"
using namespace openMVG;

#include "testing/testing.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
using namespace std;
using std::string;

static const int CARD_DESCS = 12;
static const int DESC_LENGTH = 128;
typedef Descriptor<float, DESC_LENGTH> Desc_t;

TEST(descriptorIO, ASCII) {
  // Create an input series of descriptor
  std::vector<Desc_t > vec_descs;
  for(int i = 0; i < CARD_DESCS; ++i)  {
    Desc_t desc;
    for (int j = 0; j < DESC_LENGTH; ++j)
      desc[j] = i*DESC_LENGTH+j;
    vec_descs.push_back(desc);
  }

  //Save them to a file
  saveDescsToFile("tempDescs.desc", vec_descs);

  //Read the saved data and compare to input (to check write/read IO)
  std::vector<Desc_t > vec_descs_read;
  loadDescsFromFile("tempDescs.desc", vec_descs_read);
  EXPECT_EQ(CARD_DESCS, vec_descs_read.size());

  for(int i = 0; i < CARD_DESCS; ++i) {
    for (int j = 0; j < DESC_LENGTH; ++j)
      EXPECT_EQ(vec_descs[i][j], vec_descs_read[i][j]);
  }
}

//Test binary export of descriptor
TEST(descriptorIO, BINARY) {
  // Create an input series of descriptor
  std::vector<Desc_t > vec_descs;
  for(int i = 0; i < CARD_DESCS; ++i)
  {
    Desc_t desc;
    for (int j = 0; j < DESC_LENGTH; ++j)
      desc[j] = i*DESC_LENGTH+j;
    vec_descs.push_back(desc);
  }

  //Save them to a file
  saveDescsToBinFile("tempDescsBin.desc", vec_descs);

  //Read the saved data and compare to input (to check write/read IO)
  std::vector<Desc_t > vec_descs_read;
  loadDescsFromBinFile("tempDescsBin.desc", vec_descs_read);
  EXPECT_EQ(CARD_DESCS, vec_descs_read.size());

  for(int i = 0; i < CARD_DESCS; ++i) {
    for (int j = 0; j < DESC_LENGTH; ++j)
      EXPECT_EQ(vec_descs[i][j], vec_descs_read[i][j]);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
