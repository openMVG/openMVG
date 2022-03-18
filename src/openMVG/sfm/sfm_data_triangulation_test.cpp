// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//-----------------
// Test summary:
//-----------------
// - Test sfm_data track triangulation
//-----------------

#include "openMVG/sfm/pipelines/pipelines_test.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include "testing/testing.h"

using namespace openMVG;
using namespace openMVG::sfm;

TEST(SFM_DATA_TRIANGULATION, BLIND) {

  const int nviews = 6;
  const int npoints = 32;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, cameras::PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  // Set fake data for landmark position
  for (auto& obs: sfm_data_2.structure)
  {
    obs.second.X.fill(0);
  }

  // Test the blind triangulation engine
  SfM_Data_Structure_Computation_Blind triangulation_engine;

  // Triangulate the landmarks and look if the 3D point are recovered
  {
    triangulation_engine.triangulate(sfm_data_2);
    EXPECT_EQ(npoints, sfm_data_2.structure.size());
    for (const auto& obs: sfm_data_2.structure)
    {
      EXPECT_MATRIX_NEAR(sfm_data.structure.at(obs.first).X, obs.second.X, 1e-8);
    }
  }

  // A landmark with empty information should be removed by the triangulation engine
  sfm_data_2.structure[0].obs.clear();
  {
    triangulation_engine.triangulate(sfm_data_2);
    EXPECT_EQ(npoints - 1, sfm_data_2.structure.size());
  }
}

TEST(SFM_DATA_TRIANGULATION, ROBUST) {

  const int nviews = 6;
  const int npoints = 32;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  // Translate the input dataset to a SfM_Data scene
  const SfM_Data sfm_data = getInputScene(d, config, cameras::PINHOLE_CAMERA);

  // Remove poses and structure
  SfM_Data sfm_data_2 = sfm_data;
  // Set fake data for landmark position
  for (auto& obs: sfm_data_2.structure)
  {
    obs.second.X.fill(0);
  }

  // Test the robust triangulation engine
  SfM_Data_Structure_Computation_Robust triangulation_engine;

  // Triangulate the landmarks and look if the 3D point are recovered
  {
    triangulation_engine.triangulate(sfm_data_2);
    EXPECT_EQ(npoints, sfm_data_2.structure.size());
    for (const auto& obs: sfm_data_2.structure)
    {
      EXPECT_MATRIX_NEAR(sfm_data.structure.at(obs.first).X, obs.second.X, 1e-8);
    }
  }


  // A landmark with empty information should be removed by the triangulation engine
  sfm_data_2.structure[0].obs.clear();
  {
    triangulation_engine.triangulate(sfm_data_2);
    EXPECT_EQ(npoints - 1, sfm_data_2.structure.size());
  }

  // A landmark with an outlier should still be there if triangulable
  sfm_data_2.structure.clear();
  sfm_data_2.structure[0] = sfm_data.structure.at(0);
  {
    sfm_data_2.structure[0].obs[0].x += Vec2(-10, -10000);
    sfm_data_2.structure[0].obs[1].x += Vec2(-10000, -10000);
    triangulation_engine.triangulate(sfm_data_2);
    EXPECT_EQ(1, sfm_data_2.structure.size());
  }

}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
