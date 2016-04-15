// Copyright (c) 2013,2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/half_space_intersection.hpp"
#include "openMVG/geometry/frustum.hpp"
#include "openMVG/geometry/box.hpp"

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/projection.hpp"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>

using namespace openMVG;
using namespace openMVG::geometry;
using namespace openMVG::geometry::halfPlane;
using namespace std;

// Export defined frustum in PLY file for viewing
bool export_Ply(const Frustum & frustum, const std::string & filename);

// Export defined Box in PLY file for viewing
bool export_Ply(const Box & box, const std::string & filename, const Vec3 color = Vec3(0, 255, 0));

//--
// Box/Camera frustum intersection unit test
//--

TEST(box_frustum, intersection)
{
  const int focal = 1000;
  const int principal_Point = 500;
  //-- Setup a circular camera rig or "cardioid".
  const int iNviews = 4;
  const int iNbPoints = 6;
  const NViewDataSet d =
    NRealisticCamerasRing(
    iNviews, iNbPoints,
    nViewDatasetConfigurator(focal, focal, principal_Point, principal_Point, 5, 0));


  const Box box(Vec3::Zero(), sqrt(1.));
  {
    std::ostringstream os;
    os << "box.ply";
    export_Ply(box, os.str());
  }

  // Test with infinite Frustum for each camera
  {
    std::vector<Frustum> vec_frustum;
    for (int i=0; i < iNviews; ++i)
    {
      const Frustum f (principal_Point*2, principal_Point*2, d._K[i], d._R[i], d._C[i]);
      EXPECT_TRUE(f.intersect(box));
      EXPECT_TRUE(box.intersect(f));

      std::ostringstream os;
      os << i << "frust.ply";
      export_Ply(f, os.str());
    }
  }

  // Test with truncated frustum
  {
    // Build frustum with near and far plane defined by min/max depth per camera
    for (int i=0; i < iNviews; ++i)
    {
      double minDepth = std::numeric_limits<double>::max();
      double maxDepth = std::numeric_limits<double>::min();
      for (int j=0; j < iNbPoints; ++j)
      {
        const double depth = Depth(d._R[i], d._t[i], d._X.col(j));
        if (depth < minDepth)
          minDepth = depth;
        if (depth > maxDepth)
          maxDepth = depth;
      }
      const Frustum f(principal_Point*2, principal_Point*2,
          d._K[i], d._R[i], d._C[i], minDepth, maxDepth);
      
      EXPECT_TRUE(f.intersect(box));
      EXPECT_TRUE(box.intersect(f));
    }    
  }
}

TEST(box_frustum, no_intersection)
{
  const int focal = 1000;
  const int principal_Point = 500;
  //-- Setup a circular camera rig or "cardioid".
  const int iNviews = 4;
  const int iNbPoints = 6;
  const NViewDataSet d =
    NRealisticCamerasRing(
      iNviews, iNbPoints,
      nViewDatasetConfigurator(focal, focal, principal_Point, principal_Point, 5, 0));

  // Put the box out of field of the camera
  // (since camera are Y up, we move the box along Y axis)
  const Vec3 position(0, -4, 0);
  const Box box(position, sqrt(1.));
  {
    std::ostringstream os;
    os << "box.ply";
    export_Ply(box, os.str());
  }

  // Test with infinite Frustum for each camera
  {
    std::vector<Frustum> vec_frustum;
    for (int i = 0; i < iNviews; ++i)
    {
      const Frustum f(principal_Point * 2, principal_Point * 2, d._K[i], d._R[i], d._C[i]);
      EXPECT_FALSE(f.intersect(box));
      EXPECT_FALSE(box.intersect(f));

      std::ostringstream os;
      os << i << "frust.ply";
      export_Ply(f, os.str());
    }
  }

  // Test with truncated frustum
  {
    // Build frustum with near and far plane defined by min/max depth per camera
    for (int i = 0; i < iNviews; ++i)
    {
      double minDepth = std::numeric_limits<double>::max();
      double maxDepth = std::numeric_limits<double>::min();
      for (int j = 0; j < iNbPoints; ++j)
      {
        const double depth = Depth(d._R[i], d._t[i], d._X.col(j));
        if (depth < minDepth)
          minDepth = depth;
        if (depth > maxDepth)
          maxDepth = depth;
      }
      const Frustum f(principal_Point * 2, principal_Point * 2,
        d._K[i], d._R[i], d._C[i], minDepth, maxDepth);

      EXPECT_FALSE(f.intersect(box));
      EXPECT_FALSE(box.intersect(f));
    }
  }
}

// Export defined frustum in PLY file for viewing
bool export_Ply(const Frustum & frustum, const std::string & filename)
{
  std::ofstream of(filename.c_str());
  if (!of.is_open())
    return false;
  // Vertex count evaluation
  // Faces count evaluation
  const size_t vertex_count = frustum.frustum_points().size();
  const size_t face_count = frustum.isInfinite() ? 4 : 6;

  of << "ply" << '\n'
    << "format ascii 1.0" << '\n'
    << "element vertex " << vertex_count << '\n'
    << "property float x" << '\n'
    << "property float y" << '\n'
    << "property float z" << '\n'
    << "element face " << face_count << '\n'
    << "property list uchar int vertex_index" << '\n'
    << "end_header" << '\n';

  // Export frustums points
  {
    const std::vector<Vec3> & points = frustum.frustum_points();
    for (int i = 0; i < points.size(); ++i)
      of << points[i].transpose() << '\n';
  }

  // Export frustums faces
  {
    if (frustum.isInfinite())
    {
      // infinite frustum: drawn normalized cone: 4 faces
      of
        << "3 " << 0 << ' ' << 4 << ' ' << 1 << '\n'
        << "3 " << 0 << ' ' << 1 << ' ' << 2 << '\n'
        << "3 " << 0 << ' ' << 2 << ' ' << 3 << '\n'
        << "3 " << 0 << ' ' << 3 << ' ' << 4 << '\n';
    }
    else
    {
      of
        << "4 " << 0 << ' ' << 1 << ' ' << 2 << ' ' << 3 << '\n'
        << "4 " << 0 << ' ' << 1 << ' ' << 5 << ' ' << 4 << '\n'
        << "4 " << 1 << ' ' << 5 << ' ' << 6 << ' ' << 2 << '\n'
        << "4 " << 3 << ' ' << 7 << ' ' << 6 << ' ' << 2 << '\n'
        << "4 " << 0 << ' ' << 4 << ' ' << 7 << ' ' << 3 << '\n'
        << "4 " << 4 << ' ' << 5 << ' ' << 6 << ' ' << 7 << '\n';
    }
  }
  of.flush();
  const bool bOk = of.good();
  of.close();
  return bOk;
}

// Export defined Box in PLY file for viewing
bool export_Ply(const Box & box, const std::string & filename, const Vec3 color)
{
  std::ofstream of(filename.c_str());
  if (!of.is_open())
    return false;
  // Vertex count evaluation
  // Faces count evaluation
  const size_t vertex_count = 8;
  const size_t face_count = 6;

  of << "ply" << '\n'
    << "format ascii 1.0" << '\n'
    << "element vertex " << vertex_count << '\n'
    << "property float x" << '\n'
    << "property float y" << '\n'
    << "property float z" << '\n'
    << "property uchar red" << '\n'
    << "property uchar green" << '\n'
    << "property uchar blue" << '\n'
    << "element face " << face_count << '\n'
    << "property list uchar int vertex_index" << '\n'
    << "end_header" << '\n';

  // Export box points
  {
    for (int i = 0; i < 8; ++i)
      of << box.points[i].transpose() << " " << color.cast<int>() << '\n';
  }

  of
    // top & bottom planes
    << "3 " << 0 << ' ' << 1 << ' ' << 3 << '\n'
    << "3 " << 4 << ' ' << 7 << ' ' << 5 << '\n'
    // remaining planes
    << "3 " << 0 << ' ' << 4 << ' ' << 1 << '\n'
    << "3 " << 0 << ' ' << 3 << ' ' << 4 << '\n'
    << "3 " << 3 << ' ' << 2 << ' ' << 7 << '\n'
    << "3 " << 1 << ' ' << 5 << ' ' << 2 << '\n';

  of.flush();
  const bool bOk = of.good();
  of.close();
  return bOk;
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
