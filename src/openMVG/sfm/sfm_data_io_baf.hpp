
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_IO_BAF_HPP
#define OPENMVG_SFM_DATA_IO_BAF_HPP

#include "openMVG/sfm/sfm_data_io.hpp"
#include <fstream>

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
static bool Save_BAF(
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
      << sfm_data.GetPoses().size() << '\n'
      << sfm_data.GetLandmarks().size() << '\n';

    const Intrinsics & intrinsics = sfm_data.GetIntrinsics();
    for (Intrinsics::const_iterator iterIntrinsic = intrinsics.begin();
      iterIntrinsic != intrinsics.end(); ++iterIntrinsic)
    {
        //get params
        const std::vector<double> intrinsicsParams = iterIntrinsic->second.get()->getParams();
        std::copy(intrinsicsParams.begin(), intrinsicsParams.end(),
          std::ostream_iterator<double>(stream, " "));
        stream << "\n";
    }

    const Poses & poses = sfm_data.GetPoses();
    for (Poses::const_iterator iterPose = poses.begin();
      iterPose != poses.end(); ++iterPose)
    {
      // [Rotation col major 3x3; camera center 3x1]
      const double * rotation = iterPose->second.rotation().data();
      std::copy(rotation, rotation+9, std::ostream_iterator<double>(stream, " "));
      const double * center = iterPose->second.center().data();
      std::copy(center, center+3, std::ostream_iterator<double>(stream, " "));
      stream << "\n";
    }

    const Landmarks & landmarks = sfm_data.GetLandmarks();
    for (Landmarks::const_iterator iterLandmarks = landmarks.begin();
      iterLandmarks != landmarks.end();
      ++iterLandmarks)
    {
      // Export visibility information
      // X Y Z #observations id_cam id_pose x y ...
      const double * X = iterLandmarks->second.X.data();
      std::copy(X, X+3, std::ostream_iterator<double>(stream, " "));
      const Observations & obs = iterLandmarks->second.obs;
      stream << obs.size() << " ";
      for (Observations::const_iterator iterOb = obs.begin();
        iterOb != obs.end(); ++iterOb)
      {
        const IndexT id_view = iterOb->first;
        const View * v = sfm_data.GetViews().at(id_view).get();
        stream << v->id_intrinsic << ' ' << v->id_pose << ' '
          << iterOb->second.x(0) << ' '
          << iterOb->second.x(1) << ' ';
      }
      stream << "\n";
    }

    stream.flush();
    bOk = stream.good();
    stream.close();
  }
  // Export View filenames as an imgList.txt file
  {
    const std::string sFile = stlplus::create_filespec(
      stlplus::folder_part(filename), "imgList", "txt");

    stream.open(sFile.c_str());
    if (!stream.is_open())
      return false;
    for (Views::const_iterator iterV = sfm_data.GetViews().begin();
      iterV != sfm_data.GetViews().end();
      ++ iterV)
    {
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        iterV->second->s_Img_path);
      stream << sView_filename << "\n";
    }
    stream.flush();
    bOk = stream.good();
    stream.close();
  }
  return bOk;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_IO_PLY_HPP
