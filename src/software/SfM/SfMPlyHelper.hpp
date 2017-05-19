// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_PLY_HELPER_H
#define OPENMVG_SFM_PLY_HELPER_H

#include "openMVG/numeric/numeric.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace openMVG{
namespace plyHelper{

/// Export 3D point vector to PLY format
inline
bool
exportToPly
(
  const std::vector<Vec3> & vec_points,
  const std::string & sFileName
)
{
  std::ofstream outfile(sFileName.c_str());
  if (!outfile.is_open())
    return false;

  outfile << "ply"
    << std::endl << "format ascii 1.0"
    << std::endl << "element vertex " << vec_points.size()
    << std::endl << "property double x"
    << std::endl << "property double y"
    << std::endl << "property double z"
    << std::endl << "property uchar red"
    << std::endl << "property uchar green"
    << std::endl << "property uchar blue"
    << std::endl << "end_header" << std::endl;

  outfile << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  for (size_t i=0; i < vec_points.size(); ++i)
  {
    outfile
      << vec_points[i](0) << ' '
      << vec_points[i](1) << ' '
      << vec_points[i](2) << ' '
      << "255 255 255" << "\n";
  }
  const bool bOk = outfile.good();
  outfile.close();
  return bOk;
}

/// Export 3D point vector and camera position to PLY format
inline bool exportToPly
(
  const std::vector<Vec3> & vec_points,
  const std::vector<Vec3> & vec_camPos,
  const std::string & sFileName,
  const std::vector<Vec3> * vec_coloredPoints = nullptr
)
{
  std::ofstream outfile(sFileName.c_str());
  if (!outfile.is_open())
    return false;

  outfile << "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << vec_points.size()+vec_camPos.size()
    << '\n' << "property double x"
    << '\n' << "property double y"
    << '\n' << "property double z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << std::endl;

  outfile << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  for (size_t i=0; i < vec_points.size(); ++i)  {
    if (vec_coloredPoints == nullptr)
      outfile
        << vec_points[i](0) << ' '
        << vec_points[i](1) << ' '
        << vec_points[i](2) << ' '
        << "255 255 255\n";
    else
      outfile
        << vec_points[i](0) << ' '
        << vec_points[i](1) << ' '
        << vec_points[i](2) << ' '
        << static_cast<int>((*vec_coloredPoints)[i](0)) << ' '
        << static_cast<int>((*vec_coloredPoints)[i](1)) << ' '
        << static_cast<int>((*vec_coloredPoints)[i](2))
        << "\n";
  }

  for (size_t i=0; i < vec_camPos.size(); ++i)  {
    outfile
      << vec_camPos[i](0) << ' '
      << vec_camPos[i](1) << ' '
      << vec_camPos[i](2) << ' '
      << "0 255 0\n";
  }
  outfile.flush();
  const bool bOk = outfile.good();
  outfile.close();
  return bOk;
}

} // namespace plyHelper
} // namespace openMVG

#endif // OPENMVG_SFM_PLY_HELPER_H

