// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/params/params.hpp"
#include "openMVG/params/params_camera_io.hpp"

#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <sstream>

using namespace openMVG;
using namespace openMVG::params;

// Create a Params file with default parameters
paramsCamera create_test_scene()
{
	paramsCamera params_data;
  return params_data;
}

TEST(Params_Camera_IO, SAVE_LOAD_JSON) {

  const std::vector<std::string> ext_Type = {"json"};

  for (int i=0; i < ext_Type.size(); ++i)
  {
    std::ostringstream os;
    os << "SAVE_LOAD" << "." << ext_Type[i];
    const std::string filename = os.str();
    std::cout << "Testing:" << filename << std::endl;

  // SAVE
  {
    const paramsCamera params_data = create_test_scene();
    EXPECT_TRUE( Save(params_data, filename) );
  }

  // LOAD
  {
    const paramsCamera params_data = create_test_scene();
    params_data.kMatrix="1538.003117;0.0;963.854798;0.0;1538.431624;455.046237;0.0;0.0;1";
    EXPECT_TRUE( Save(params_data, filename) );
    paramsCamera params_data_load;
    EXPECT_TRUE( Load(params_data_load, filename) );
    EXPECT_EQ( params_data_load.kMatrix, params_data.kMatrix);
  }

  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
