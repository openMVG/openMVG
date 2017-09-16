// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_filters_frustum.hpp"

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/stl/stl.hpp"

#include "third_party/progress/progress_display.hpp"

#include <fstream>
#include <iomanip>
#include <iterator>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::geometry::halfPlane;

// Constructor
Frustum_Filter::Frustum_Filter
(
  const SfM_Data & sfm_data,
  const double zNear,
  const double zFar,
  const NearFarPlanesT & z_near_z_far
)
: z_near_z_far_perView(z_near_z_far)
{
  //-- Init Z_Near & Z_Far for all valid views
  init_z_near_z_far_depth(sfm_data, zNear, zFar);
  const bool bComputed_Z = (zNear == -1. && zFar == -1.) && !sfm_data.structure.empty();
  _bTruncated = (zNear != -1. && zFar != -1.) || bComputed_Z;
  initFrustum(sfm_data);
}

// Init a frustum for each valid views of the SfM scene
void Frustum_Filter::initFrustum
(
  const SfM_Data & sfm_data
)
{
  for (NearFarPlanesT::const_iterator it = z_near_z_far_perView.begin();
      it != z_near_z_far_perView.end(); ++it)
  {
    const View * view = sfm_data.GetViews().at(it->first).get();
    if (!sfm_data.IsPoseAndIntrinsicDefined(view))
      continue;
    Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
    if (!isPinhole(iterIntrinsic->second->getType()))
      continue;

    const Pose3 pose = sfm_data.GetPoseOrDie(view);

    const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic->second.get());
    if (cam == nullptr)
      continue;

    if (!_bTruncated) // use infinite frustum
    {
      const Frustum f(
        cam->w(), cam->h(), cam->K(),
        pose.rotation(), pose.center());
      frustum_perView[view->id_view] = f;
    }
    else // use truncated frustum with defined Near and Far planes
    {
      const Frustum f(cam->w(), cam->h(), cam->K(),
        pose.rotation(), pose.center(), it->second.first, it->second.second);
      frustum_perView[view->id_view] = f;
    }
  }
}

Pair_Set Frustum_Filter::getFrustumIntersectionPairs
(
  const std::vector<HalfPlaneObject>& bounding_volume
)
const
{
  Pair_Set pairs;
  // List active view Id
  std::vector<IndexT> viewIds;
  viewIds.reserve(z_near_z_far_perView.size());
  std::transform(z_near_z_far_perView.begin(), z_near_z_far_perView.end(),
    std::back_inserter(viewIds), stl::RetrieveKey());

  C_Progress_display my_progress_bar(
    viewIds.size() * (viewIds.size()-1)/2,
    std::cout, "\nCompute frustum intersection\n");

  // Exhaustive comparison (use the fact that the intersect function is symmetric)
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for
#endif
  for (int i = 0; i < (int)viewIds.size(); ++i)
  {
    // Prepare vector of intersecting objects (within loop to keep it
    // thread-safe)
    std::vector<HalfPlaneObject> objects = bounding_volume;
    objects.insert(objects.end(),
                   { frustum_perView.at(viewIds[i]), HalfPlaneObject() });

    for (size_t j = i+1; j < viewIds.size(); ++j)
    {
      objects.back() = frustum_perView.at(viewIds[j]);
      if (intersect(objects))
      {
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          pairs.insert({viewIds[i], viewIds[j]});
        }
      }
      // Progress bar update
      ++my_progress_bar;
    }
  }
  return pairs;
}

// Export defined frustum in PLY file for viewing
bool Frustum_Filter::export_Ply
(
  const std::string & filename
)
const
{
  std::ofstream of(filename.c_str());
  if (!of.is_open())
    return false;
  // Vertex count evaluation
  // Faces count evaluation
  size_t vertex_count = 0;
  size_t face_count = 0;
  for (FrustumsT::const_iterator it = frustum_perView.begin();
    it != frustum_perView.end(); ++it)
  {
    if (it->second.isInfinite())
    {
      vertex_count += 5;
      face_count += 5; // 4 triangles + 1 quad
    }
    else // truncated
    {
      vertex_count += 8;
      face_count += 6; // 6 quads
    }
  }

  of << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  of << "ply" << '\n'
    << "format ascii 1.0" << '\n'
    << "element vertex " << vertex_count << '\n'
    << "property double x" << '\n'
    << "property double y" << '\n'
    << "property double z" << '\n'
    << "element face " << face_count << '\n'
    << "property list uchar int vertex_index" << '\n'
    << "end_header" << '\n';

  // Export frustums points
  for (FrustumsT::const_iterator it = frustum_perView.begin();
    it != frustum_perView.end(); ++it)
  {
    const std::vector<Vec3> & points = it->second.frustum_points();
    for (size_t i=0; i < points.size(); ++i)
      of << points[i].transpose() << '\n';
  }

  // Export frustums faces
  IndexT count = 0;
  for (FrustumsT::const_iterator it = frustum_perView.begin();
    it != frustum_perView.end(); ++it)
  {
    if (it->second.isInfinite()) // infinite frustum: drawn normalized cone: 4 faces
    {
      of << "3 " << count + 0 << ' ' << count + 1 << ' ' << count + 2 << '\n'
        << "3 " << count + 0 << ' ' << count + 2 << ' ' << count + 3 << '\n'
        << "3 " << count + 0 << ' ' << count + 3 << ' ' << count + 4 << '\n'
        << "3 " << count + 0 << ' ' << count + 4 << ' ' << count + 1 << '\n'
        << "4 " << count + 1 << ' ' << count + 2 << ' ' << count + 3 << ' ' << count + 4 << '\n';
      count += 5;
    }
    else // truncated frustum: 6 faces
    {
      of << "4 " << count + 0 << ' ' << count + 1 << ' ' << count + 2 << ' ' << count + 3 << '\n'
        << "4 " << count + 0 << ' ' << count + 1 << ' ' << count + 5 << ' ' << count + 4 << '\n'
        << "4 " << count + 1 << ' ' << count + 5 << ' ' << count + 6 << ' ' << count + 2 << '\n'
        << "4 " << count + 3 << ' ' << count + 7 << ' ' << count + 6 << ' ' << count + 2 << '\n'
        << "4 " << count + 0 << ' ' << count + 4 << ' ' << count + 7 << ' ' << count + 3 << '\n'
        << "4 " << count + 4 << ' ' << count + 5 << ' ' << count + 6 << ' ' << count + 7 << '\n';
      count += 8;
    }
  }
  of.flush();
  const bool bOk = of.good();
  of.close();
  return bOk;
}

void Frustum_Filter::init_z_near_z_far_depth
(
  const SfM_Data & sfm_data,
  const double zNear,
  const double zFar
)
{
  // If z_near & z_far are -1 and structure if not empty,
  //  compute the values for each camera and the structure
  const bool bComputed_Z = (zNear == -1. && zFar == -1.) && !sfm_data.structure.empty();
  if (bComputed_Z)  // Compute the near & far planes from the structure and view observations
  {
    for (Landmarks::const_iterator itL = sfm_data.GetLandmarks().begin();
      itL != sfm_data.GetLandmarks().end(); ++itL)
    {
      const Landmark & landmark = itL->second;
      const Vec3 & X = landmark.X;
      for (Observations::const_iterator iterO = landmark.obs.begin();
        iterO != landmark.obs.end(); ++iterO)
      {
        const IndexT id_view = iterO->first;
        const View * view = sfm_data.GetViews().at(id_view).get();
        if (!sfm_data.IsPoseAndIntrinsicDefined(view))
          continue;

        const Pose3 pose = sfm_data.GetPoseOrDie(view);
        const double z = pose.depth(X);
        NearFarPlanesT::iterator itZ = z_near_z_far_perView.find(id_view);
        if (itZ != z_near_z_far_perView.end())
        {
          if ( z < itZ->second.first)
            itZ->second.first = z;
          else
          if ( z > itZ->second.second)
            itZ->second.second = z;
        }
        else
          z_near_z_far_perView[id_view] = {z,z};
      }
    }
  }
  else
  {
    // Init the same near & far limit for all the valid views
    for (Views::const_iterator it = sfm_data.GetViews().begin();
    it != sfm_data.GetViews().end(); ++it)
    {
      const View * view = it->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;
      if (z_near_z_far_perView.find(view->id_view) == z_near_z_far_perView.end())
        z_near_z_far_perView[view->id_view] = {zNear, zFar};
    }
  }
}

} // namespace sfm
} // namespace openMVG
