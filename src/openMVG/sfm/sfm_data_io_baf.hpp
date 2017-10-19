// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_IO_BAF_HPP
#define OPENMVG_SFM_SFM_DATA_IO_BAF_HPP

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "openMVG/sfm/sfm_data_io.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG {
namespace sfm {

/// Save SfM_Data in an ASCII BAF (Bundle Adjustment File).
// --Header
// #Intrinsics
// #Poses
// #Landmarks
// --Data
// Intrinsic parameters [foc ppx ppy, ...]
// Poses [angle axis, camera center]
// Landmarks [X Y Z #observations id_intrinsic id_pose x y ...]
//--
//- Export also a _imgList.txt file with View filename and id_intrinsic & id_pose.
// filename id_intrinsic id_pose
// The ids allow to establish a link between 3D point observations & the corresponding views
//--
// Export missing poses as Identity pose to keep tracking of the original id_pose indexes
inline bool Save_BAF(
  const SfM_Data & sfm_data,
  const std::string & filename,
  ESfM_Data flags_part)
{
  std::ofstream stream(filename.c_str());
  if (!stream.is_open())
    return false;

  bool bOk = false;
  {
    stream
      << sfm_data.GetIntrinsics().size() << '\n'
      << sfm_data.GetViews().size() << '\n'
      << sfm_data.GetLandmarks().size() << '\n';

    const Intrinsics & intrinsics = sfm_data.GetIntrinsics();
    for (const auto & iterIntrinsic : intrinsics )
    {
      //get params
      const std::vector<double> intrinsicsParams = iterIntrinsic.second->getParams();
      std::copy(intrinsicsParams.begin(), intrinsicsParams.end(),
        std::ostream_iterator<double>(stream, " "));
      stream << '\n';
    }

    const Poses & poses = sfm_data.GetPoses();
    for ( const auto & iterV : sfm_data.GetViews() )
    {
      const View * view = iterV.second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
      {
        const Mat3 R = Mat3::Identity();
        const double * rotation = R.data();
        std::copy(rotation, rotation+9, std::ostream_iterator<double>(stream, " "));
        const Vec3 C = Vec3::Zero();
        const double * center = C.data();
        std::copy(center, center+3, std::ostream_iterator<double>(stream, " "));
        stream << '\n';
      }
      else
      {
        // [Rotation col major 3x3; camera center 3x1]
        const double * rotation = poses.at(view->id_pose).rotation().data();
        std::copy(rotation, rotation+9, std::ostream_iterator<double>(stream, " "));
        const double * center = poses.at(view->id_pose).center().data();
        std::copy(center, center+3, std::ostream_iterator<double>(stream, " "));
        stream << '\n';
      }
    }

    const Landmarks & landmarks = sfm_data.GetLandmarks();
    for (const auto & iterLandmarks : landmarks )
    {
      // Export visibility information
      // X Y Z #observations id_cam id_pose x y ...
      const double * X = iterLandmarks.second.X.data();
      std::copy(X, X+3, std::ostream_iterator<double>(stream, " "));
      const Observations & obs = iterLandmarks.second.obs;
      stream << obs.size() << " ";
      for ( const auto & iterOb : obs )
      {
        const IndexT id_view = iterOb.first;
        const View * v = sfm_data.GetViews().at(id_view).get();
        stream
          << v->id_intrinsic << ' '
          << v->id_pose << ' '
          << iterOb.second.x(0) << ' ' << iterOb.second.x(1) << ' ';
      }
      stream << '\n';
    }

    stream.flush();
    bOk = stream.good();
    stream.close();
  }

  // Export View filenames & ids as an imgList.txt file
  {
    const std::string sFile = stlplus::create_filespec(
      stlplus::folder_part(filename), stlplus::basename_part(filename) + std::string("_imgList"), "txt");

    stream.open(sFile.c_str());
    if (!stream.is_open())
      return false;
    for ( const auto & iterV : sfm_data.GetViews() )
    {
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        iterV.second->s_Img_path);
      stream
        << sView_filename
        << ' ' << iterV.second->id_intrinsic
        << ' ' << iterV.second->id_pose << "\n";
    }
    stream.flush();
    bOk = stream.good();
    stream.close();
  }
  return bOk;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_IO_BAF_HPP
