// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2026 Pierre MOULON, Romuald PERROT.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_report.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/pipelines/pipelines_test.hpp"
#include "openMVG/system/logger.hpp"

#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <string>
#include <fstream>
#include <sstream>

using namespace openMVG;
using namespace openMVG::sfm;

TEST(SFM_REPORT, EmptyScene)
{
  SfM_Data empty_data;
  const std::string report_file = "empty_report.html";
  stlplus::file_delete(report_file);
  EXPECT_FALSE(Generate_SfM_Report(empty_data, report_file));
  EXPECT_FALSE(stlplus::file_exists(report_file));
}

TEST(SFM_REPORT, SyntheticScene)
{
  const int nviews = 6;
  const int npoints = 32;
  const nViewDatasetConfigurator config;
  const NViewDataSet d = NRealisticCamerasRing(nviews, npoints, config);

  const SfM_Data sfm_data = getInputScene(d, config, cameras::PINHOLE_CAMERA);

  const std::string report_file = "synthetic_report.html";
  stlplus::file_delete(report_file);
  EXPECT_TRUE(Generate_SfM_Report(sfm_data, report_file));
  EXPECT_TRUE(stlplus::file_exists(report_file));

  // Check generated SVG files
  const std::string folder = stlplus::folder_part(report_file);
  const std::string tracks_svg = stlplus::create_filespec(folder, "tracks_histogram", "svg");
  const std::string residuals_svg = stlplus::create_filespec(folder, "residuals_histogram", "svg");
  const std::string parallax_svg = stlplus::create_filespec(folder, "parallax_angles_histogram", "svg");

  EXPECT_TRUE(stlplus::file_exists(tracks_svg));
  EXPECT_TRUE(stlplus::file_exists(residuals_svg));
  EXPECT_TRUE(stlplus::file_exists(parallax_svg));

  // Read HTML content and verify key statistics tables are present
  std::ifstream f(report_file);
  EXPECT_TRUE(f.is_open());
  std::stringstream buffer;
  buffer << f.rdbuf();
  const std::string html = buffer.str();

  EXPECT_TRUE(html.find("Track length statistics:") != std::string::npos);
  EXPECT_TRUE(html.find("Track Length Categories") != std::string::npos);
  EXPECT_TRUE(html.find("Detailed Distribution per Length") != std::string::npos);
  EXPECT_TRUE(html.find("Triangulation Parallax Angles") != std::string::npos);
  EXPECT_TRUE(html.find("Camera Covisibility Graph:") != std::string::npos);
  EXPECT_TRUE(html.find("Reprojection Error (2D Norm) Statistics:") != std::string::npos);

  // Clean up
  stlplus::file_delete(report_file);
  stlplus::file_delete(tracks_svg);
  stlplus::file_delete(residuals_svg);
  stlplus::file_delete(parallax_svg);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
/* ************************************************************************* */
