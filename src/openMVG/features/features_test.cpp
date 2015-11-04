// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/features.hpp"
using namespace openMVG;
using namespace openMVG::features;

#include "testing/testing.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
using namespace std;
using std::string;

//-- Test features
static const int CARD = 12;

typedef SIOPointFeature Feature_T;
typedef std::vector<Feature_T> Feats_T;

TEST(featureIO, ASCII) {
  Feats_T vec_feats;
  for(int i = 0; i < CARD; ++i)  {
    vec_feats.push_back(Feature_T(i, i*2, i*3, i*4));
  }

  //Save them to a file
  saveFeatsToFile("tempFeats.feat", vec_feats);

  //Read the saved data and compare to input (to check write/read IO)
  Feats_T vec_feats_read;
  loadFeatsFromFile("tempFeats.feat", vec_feats_read);
  EXPECT_EQ(CARD, vec_feats_read.size());

  for(int i = 0; i < CARD; ++i) {
    EXPECT_EQ(vec_feats[i], vec_feats_read[i]);
    EXPECT_EQ(vec_feats[i].coords(), vec_feats_read[i].coords());
    EXPECT_EQ(vec_feats[i].scale(), vec_feats_read[i].scale());
    EXPECT_EQ(vec_feats[i].orientation(), vec_feats_read[i].orientation());
  }
}

//-- Test descriptors

static const int DESC_LENGTH = 128;
typedef Descriptor<float, DESC_LENGTH> Desc_T;
typedef std::vector<Desc_T> Descs_T;

TEST(descriptorIO, ASCII) {
  // Create an input series of descriptor
  Descs_T vec_descs;
  for(int i = 0; i < CARD; ++i)  {
    Desc_T desc;
    for (int j = 0; j < DESC_LENGTH; ++j)
      desc[j] = i*DESC_LENGTH+j;
    vec_descs.push_back(desc);
  }

  //Save them to a file
  saveDescsToFile("tempDescs.desc", vec_descs);

  //Read the saved data and compare to input (to check write/read IO)
  Descs_T vec_descs_read;
  loadDescsFromFile("tempDescs.desc", vec_descs_read);
  EXPECT_EQ(CARD, vec_descs_read.size());

  for(int i = 0; i < CARD; ++i) {
    for (int j = 0; j < DESC_LENGTH; ++j)
      EXPECT_EQ(vec_descs[i][j], vec_descs_read[i][j]);
  }
}

//Test binary export of descriptor
TEST(descriptorIO, BINARY) {
  // Create an input series of descriptor
  Descs_T vec_descs;
  for(int i = 0; i < CARD; ++i)
  {
    Desc_T desc;
    for (int j = 0; j < DESC_LENGTH; ++j)
      desc[j] = i*DESC_LENGTH+j;
    vec_descs.push_back(desc);
  }

  //Save them to a file
  saveDescsToBinFile("tempDescsBin.desc", vec_descs);

  //Read the saved data and compare to input (to check write/read IO)
  Descs_T vec_descs_read;
  loadDescsFromBinFile("tempDescsBin.desc", vec_descs_read);
  EXPECT_EQ(CARD, vec_descs_read.size());

  for(int i = 0; i < CARD; ++i) {
    for (int j = 0; j < DESC_LENGTH; ++j)
      EXPECT_EQ(vec_descs[i][j], vec_descs_read[i][j]);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
