
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <sstream>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;

// Test a border case (empty scene)
TEST(SfM_Data_IntrinsicGrouping, Empty)
{
  SfM_Data sfm_data;
  GroupSharedIntrinsics(sfm_data);
}

// Two similar intrinsics object must be grouped as one
TEST(SfM_Data_IntrinsicGrouping, Grouping_One)
{
  SfM_Data sfm_data;
  GroupSharedIntrinsics(sfm_data);
  // One view, one intrinsic
  sfm_data.views[0] = std::make_shared<View>("", 0, 0, 0);
  sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic>(0);
  // One view, one intrinsic
  sfm_data.views[1] = std::make_shared<View>("", 1, 1, 1);
  sfm_data.intrinsics[1] = std::make_shared<Pinhole_Intrinsic>(0);
  
  CHECK_EQUAL(2, sfm_data.intrinsics.size());
  GroupSharedIntrinsics(sfm_data);
  CHECK_EQUAL(1, sfm_data.intrinsics.size());
  CHECK_EQUAL(0, sfm_data.views[0]->id_intrinsic);
  CHECK_EQUAL(0, sfm_data.views[1]->id_intrinsic); 
}

// Unique intrinsics object must not be grouped
TEST(SfM_Data_IntrinsicGrouping, No_Grouping)
{
  SfM_Data sfm_data;
  GroupSharedIntrinsics(sfm_data);
  const int nbView = 10;
  for(int i=0; i < nbView; ++i)
  {
    // Add one view, one intrinsic
    sfm_data.views[i] = std::make_shared<View>("", i, i, i);
    sfm_data.intrinsics[i] = std::make_shared<Pinhole_Intrinsic>(i);
  }
  
  CHECK_EQUAL(nbView, sfm_data.intrinsics.size());
  GroupSharedIntrinsics(sfm_data);
  CHECK_EQUAL(nbView, sfm_data.intrinsics.size());
}

// Similar intrinsics object must be grouped as one
TEST(SfM_Data_IntrinsicGrouping, Grouping_Two)
{
  SfM_Data sfm_data;
  GroupSharedIntrinsics(sfm_data);
  // Define separate intrinsics that share common properties
  // first block of intrinsics
  {
    sfm_data.views[0] = std::make_shared<View>("", 0, 0, 0);
    sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic>(0);
    sfm_data.views[1] = std::make_shared<View>("", 1, 1, 1);
    sfm_data.intrinsics[1] = std::make_shared<Pinhole_Intrinsic>(0);
  }
  // second block of intrinsics
  {
    sfm_data.views[2] = std::make_shared<View>("", 2, 2, 2);
    sfm_data.intrinsics[2] = std::make_shared<Pinhole_Intrinsic>(1);
    sfm_data.views[3] = std::make_shared<View>("", 3, 3, 3);
    sfm_data.intrinsics[3] = std::make_shared<Pinhole_Intrinsic>(1);
  }
  
  CHECK_EQUAL(4, sfm_data.intrinsics.size());
  GroupSharedIntrinsics(sfm_data);
  CHECK_EQUAL(2, sfm_data.intrinsics.size());
  CHECK_EQUAL(0, sfm_data.views[0]->id_intrinsic);
  CHECK_EQUAL(0, sfm_data.views[1]->id_intrinsic);
  CHECK_EQUAL(1, sfm_data.views[2]->id_intrinsic);
  CHECK_EQUAL(1, sfm_data.views[3]->id_intrinsic);
}

// Similar intrinsics object must be grouped as one (test with various camera type)
TEST(SfM_Data_IntrinsicGrouping, Grouping_Two_Different_Camera_Type)
{
  SfM_Data sfm_data;
  GroupSharedIntrinsics(sfm_data);
  // Define separate intrinsics that share common properties
  // first block of intrinsics
  {
    sfm_data.views[0] = std::make_shared<View>("", 0, 0, 0);
    sfm_data.intrinsics[0] = std::make_shared<Pinhole_Intrinsic>(0);
    sfm_data.views[1] = std::make_shared<View>("", 1, 1, 1);
    sfm_data.intrinsics[1] = std::make_shared<Pinhole_Intrinsic>(0);
  }
  // second block of intrinsics (different type)
  {
    sfm_data.views[2] = std::make_shared<View>("", 2, 2, 2);
    sfm_data.intrinsics[2] = std::make_shared<Pinhole_Intrinsic_Radial_K1>(0);
    sfm_data.views[3] = std::make_shared<View>("", 3, 3, 3);
    sfm_data.intrinsics[3] = std::make_shared<Pinhole_Intrinsic_Radial_K1>(0);
  }
  
  CHECK_EQUAL(4, sfm_data.intrinsics.size());
  GroupSharedIntrinsics(sfm_data);
  CHECK_EQUAL(2, sfm_data.intrinsics.size());
  CHECK_EQUAL(0, sfm_data.views[0]->id_intrinsic);
  CHECK_EQUAL(0, sfm_data.views[1]->id_intrinsic);
  CHECK_EQUAL(1, sfm_data.views[2]->id_intrinsic);
  CHECK_EQUAL(1, sfm_data.views[3]->id_intrinsic);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
