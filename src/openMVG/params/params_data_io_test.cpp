// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/params/params.hpp"
#include "openMVG/params/params_data_io.hpp"

#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <sstream>

using namespace openMVG;
using namespace openMVG::params;

// Create a Params file with default parameters
paramsData create_test_scene()
{
  paramsData params_data;
  return params_data;
}

TEST(Params_Data_IO, SAVE_LOAD_JSON) {

  const std::vector<std::string> ext_Type = {"json"};

  for (int i=0; i < ext_Type.size(); ++i)
  {
    std::ostringstream os;
    os << "SAVE_LOAD" << "." << ext_Type[i];
    const std::string filename = os.str();
    std::cout << "Testing:" << filename << std::endl;

  // SAVE
  {
    const paramsData params_data = create_test_scene();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsData params_data = create_test_scene();
    EXPECT_TRUE( Save(params_data, filename) );
    paramsData params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.matching.max_matching_dist_ratio, params_data.matching.max_matching_dist_ratio);
  }

  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
