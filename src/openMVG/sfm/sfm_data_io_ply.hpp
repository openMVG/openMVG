
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_IO_PLY_HPP
#define OPENMVG_SFM_DATA_IO_PLY_HPP

#include "openMVG/sfm/sfm_data_io.hpp"
#include <fstream>

namespace openMVG
{

/// Save the structure and camera positions of a SfM_Data container as 3D points in a PLY ASCII file.
static bool Save_PLY(
  const SfM_Data & data,
  const std::string & filename,
  ESfM_Data flags_part)
{
  const bool b_structure = (flags_part & STRUCTURE) == STRUCTURE;
  const bool b_extrinsics = (flags_part & EXTRINSICS) == EXTRINSICS;

  if (!(b_structure || b_extrinsics))
    return false;

  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str());
  if (!stream.is_open())
    return false;

  bool bOk = false;
  {
    stream << "ply"
      << '\n' << "format ascii 1.0"
      << '\n' << "element vertex "
        << ((b_structure ? data.getLandmarks().size() : 0) +
           (b_extrinsics ? data.getPoses().size() : 0))
      << '\n' << "property float x"
      << '\n' << "property float y"
      << '\n' << "property float z"
      << '\n' << "property uchar red"
      << '\n' << "property uchar green"
      << '\n' << "property uchar blue"
      << '\n' << "end_header" << std::endl;

      if (b_extrinsics)
      {
        const Poses & poses = data.getPoses();
        for (Poses::const_iterator iterPose = poses.begin();
          iterPose != poses.end(); ++iterPose)  {
            stream << iterPose->second.center().transpose()
              << " 0 255 0" << "\n";
        }
      }

      if (b_structure)
      {
        const Landmarks & landmarks = data.getLandmarks();
        for (Landmarks::const_iterator iterLandmarks = landmarks.begin();
          iterLandmarks != landmarks.end();
          ++iterLandmarks)  {
          stream << iterLandmarks->second.X.transpose() << " 255 255 255" << "\n";
        }
      }
      stream.flush();
      bOk = stream.good();
      stream.close();
  }
  return bOk;
}

} // namespace openMVG

#endif // OPENMVG_SFM_DATA_IO_PLY_HPP
