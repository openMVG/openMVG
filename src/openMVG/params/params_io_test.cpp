// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/params/params_io.hpp>
#include "openMVG/params/params.hpp"
#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <sstream>

using namespace openMVG;
using namespace openMVG::params;

// Create a Params file with default parameters
paramsCamera create_test_camera()
{
  paramsCamera params_data;
  return params_data;
}
paramsDetection create_test_detection()
{
  paramsDetection params_data;
  return params_data;
}
paramsMatching create_test_matching()
{
  paramsMatching params_data;
  return params_data;
}
paramsIncrementalSfM create_test_incrementalSfM()
{
  paramsIncrementalSfM params_data;
  return params_data;
}
paramsGlobalSfM create_test_globalSfM()
{
  paramsGlobalSfM params_data;
  return params_data;
}

TEST(Params_IO, SAVE_LOAD_JSON) {

  const std::vector<std::string> ext_Type = {"json"};

  for (int i=0; i < ext_Type.size(); ++i)
  {
    std::ostringstream os;
    os << "SAVE_LOAD" << "." << ext_Type[i];
    const std::string filename = os.str();
    std::cout << "Testing:" << filename << std::endl;

  // Camera check
  // SAVE
  {
    const paramsCamera params_data = create_test_camera();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsCamera params_data = create_test_camera();
    EXPECT_TRUE( Save(params_data, filename) );
    paramsCamera params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.kMatrix, params_data.kMatrix);
  }

  // Detection check
  // SAVE
  {
    const paramsDetection params_data = create_test_detection();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsDetection params_data = create_test_detection();
    EXPECT_TRUE( Save(params_data, filename) );
    paramsDetection params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.feature_preset, params_data.feature_preset);
  }

  // Matching check
  // SAVE
  {
    const paramsMatching params_data = create_test_matching();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsMatching params_data = create_test_matching();
    EXPECT_TRUE( Save(params_data, filename) );
    paramsMatching params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.max_matching_dist_ratio, params_data.max_matching_dist_ratio);
  }

  // IncrementalSfM check
  // SAVE
  {
    const paramsIncrementalSfM params_data = create_test_incrementalSfM();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsIncrementalSfM params_data = create_test_incrementalSfM();
    EXPECT_TRUE( Save(params_data, filename) );
    paramsIncrementalSfM params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.init_pair_best_of_k, params_data.init_pair_best_of_k);
  }

  // GlobalSfM check
  // SAVE
  {
    const paramsGlobalSfM params_data = create_test_globalSfM();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsGlobalSfM params_data = create_test_globalSfM();
    EXPECT_TRUE( Save(params_data, filename) );
    paramsGlobalSfM params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.rotationAveragingMethod, params_data.rotationAveragingMethod);
  }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
