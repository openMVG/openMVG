// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_report.hpp"
#include "testing/testing.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <iterator>
#include <map>
#include <string>

using namespace openMVG;
using namespace openMVG::sfm;

namespace {

void AddView(SfM_Data & scene, IndexT id, const Vec3 & center)
{
  scene.views[id] = std::make_shared<View>("image.jpg", id, 0, id, 1600, 1200);
  scene.poses[id] = geometry::Pose3(Mat3::Identity(), center);
}

void AddTrack(SfM_Data & scene, const std::vector<IndexT> & views)
{
  Landmark track;
  track.X = Vec3(0, 0, 10);
  for (const IndexT id : views)
  {
    track.obs[id] = Observation(
      scene.intrinsics.at(0)->project(scene.poses.at(id)(track.X)), 0);
  }
  scene.structure[static_cast<IndexT>(scene.structure.size())] = std::move(track);
}

std::string ReadReport(const std::string & filename)
{
  std::ifstream file(filename);
  return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
}

void RemoveReport(const std::string & filename)
{
  stlplus::file_delete(filename);
  const std::string folder = stlplus::folder_part(filename);
  for (const std::string & name : {"tracks_histogram", "residuals_histogram", "parallax_angles_histogram"})
    stlplus::file_delete(stlplus::create_filespec(folder, name, "svg"));
  stlplus::folder_delete(folder);
}

} // namespace

TEST(SFM_REPORT_STATISTICS, LongAndShortTracksKeepPairStatistics)
{
  SfM_Data scene;
  scene.intrinsics[0] = std::make_shared<cameras::Pinhole_Intrinsic>(1600, 1200, 900, 800, 600);
  AddView(scene, 7, Vec3(-10, 0, 0));
  AddView(scene, 19, Vec3(10, 0, 0));
  AddView(scene, 100, Vec3(-10, 0, 0));
  AddView(scene, 1000000000, Vec3(0, 0, 0));
  for (int i = 0; i < 30; ++i)
    AddTrack(scene, {7, 19, 100, 1000000000});
  AddTrack(scene, {7, 19});
  AddTrack(scene, {7, 100}); // Coincident camera centers give zero parallax.
  AddTrack(scene, {19}); // A single observation contributes no pair or angle.

  // Observations with missing views or missing poses must not add valid pairs.
  scene.views[888] = std::make_shared<View>("unposed.jpg", 888, 0, UndefinedIndexT, 1600, 1200);
  scene.structure.at(0).obs[888] = Observation(Vec2::Zero(), 0);
  scene.structure.at(0).obs[999] = Observation(Vec2::Zero(), 0);

  const std::string folder = "sfm_report_statistics_pairs";
  stlplus::folder_create(folder);
  const std::string filename = stlplus::create_filespec(folder, "report", "html");
  EXPECT_TRUE(Generate_SfM_Report(scene, filename));
  const std::string html = ReadReport(filename);
  EXPECT_TRUE(html.find("Covisible camera pairs (sharing >= 1 point)</td><td>6</td>") != std::string::npos);
  EXPECT_TRUE(html.find("Covisible camera pairs sharing >= 30 points</td><td>6</td>") != std::string::npos);
  EXPECT_TRUE(html.find("<td>30 / 30.00 / 30.33 / 31</td>") != std::string::npos);
  // Thirty-one tracks have 90-degree maximum parallax; one has zero.
  EXPECT_TRUE(html.find("<td>0.00&deg;</td><td>87.19&deg;</td><td>90.00&deg;</td>") != std::string::npos);

  for (auto & entry : scene.structure)
  {
    const std::map<IndexT, Observation> previous(entry.second.obs.begin(), entry.second.obs.end());
    entry.second.obs.clear();
    for (auto observation = previous.rbegin(); observation != previous.rend(); ++observation)
      entry.second.obs.emplace(*observation);
  }
  EXPECT_TRUE(Generate_SfM_Report(scene, filename));
  EXPECT_EQ(html, ReadReport(filename));
  RemoveReport(filename);
}

TEST(SFM_REPORT_STATISTICS, OppositeRaysKeepMaximumParallax)
{
  SfM_Data scene;
  scene.intrinsics[0] = std::make_shared<cameras::Pinhole_Intrinsic>(1600, 1200, 900, 800, 600);
  AddView(scene, 2, Vec3(0, 0, 0));
  AddView(scene, 11, Vec3(0, 0, 20));
  AddView(scene, 28, Vec3(10, 0, 0));
  AddTrack(scene, {2, 11, 28});

  const std::string folder = "sfm_report_statistics_opposite";
  stlplus::folder_create(folder);
  const std::string filename = stlplus::create_filespec(folder, "report", "html");
  EXPECT_TRUE(Generate_SfM_Report(scene, filename));
  const std::string html = ReadReport(filename);
  EXPECT_TRUE(html.find("<td>180.00&deg;</td><td>180.00&deg;</td><td>180.00&deg;</td>") != std::string::npos);
  EXPECT_TRUE(html.find("Covisible camera pairs (sharing >= 1 point)</td><td>3</td>") != std::string::npos);
  EXPECT_TRUE(html.find("<td>1 / 1.00 / 1.00 / 1</td>") != std::string::npos);
  RemoveReport(filename);
}

int main() { TestResult result; return TestRegistry::runAllTests(result); }
