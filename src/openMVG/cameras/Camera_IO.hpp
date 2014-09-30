
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_IO_H
#define OPENMVG_CAMERA_IO_H

#include "openMVG/cameras/PinholeCamera.hpp"
#include "openMVG/cameras/BrownPinholeCamera.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>

namespace openMVG{

/// Save a Pinhole camera to a file as a P matrix (12 doubles as binary values)
static bool save(
  const std::string & scameraFile,
  const PinholeCamera & cam)
{
  if (stlplus::extension_part(scameraFile) != "bin")
    return false;

  const Mat34 & PMat = cam._P;
  std::ofstream file( scameraFile.c_str(), std::ios::out|std::ios::binary);
  file.write((const char*)PMat.data(),(std::streamsize)(3*4)*sizeof(double));

  bool bOk = (!file.fail());
  file.close();
  return bOk;
}

/// Load a Pinhole camera from a file (read 12 doubles saved in binary values)
static bool load(
  const std::string & scameraFile,
  PinholeCamera & cam)
{
  std::vector<double> val;

  if (stlplus::extension_part(scameraFile) == "bin")
  {
    std::ifstream in(scameraFile.c_str(), std::ios::in|std::ios::binary);
    if (!in.is_open())	{
      std::cerr << "Error: failed to open file '" << scameraFile
        << "' for reading" << std::endl;
      return false;
    }
    val.resize(12);
    in.read((char*)&val[0],(std::streamsize)12*sizeof(double));
    if (in.fail())  {
      val.clear();
    }
  }
  else
    return false;

  if (val.size() == 3*4) //P Matrix
  {
    Mat34 P;
    P << val[0], val[3], val[6], val[9],
      val[1], val[4], val[7], val[10],
      val[2], val[5], val[8], val[11];

    Mat3 R,K;
    Vec3 t;
    KRt_From_P(P, &K, &R, &t);
    cam = PinholeCamera(K, R, t);
    return true;
  }
  return false;
}

/// Load a BrownPinholeCamera saved as ascii file
static bool load(
  const std::string & scameraFile,
  BrownPinholeCamera & cam)
{
  if (stlplus::extension_part(scameraFile) != "txt")
    return false;

  std::ifstream in(scameraFile.c_str(),  std::ios::in);
  if (!in.is_open())  {
    std::cerr << "Error: failed to open file '" << scameraFile
      << "' for reading" << std::endl;
    return false;
  }
  double f, ppx, ppy, k1, k2, k3;
  Mat3 R;
  Vec3 t;
  in >> f >> ppx >> ppy >> k1 >> k2 >> k3;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      in >> R(i,j);
    }
  }

  for (int i=0; i<3; ++i) {
    in >> t(i);
  }
  in.close();

  if (!in.good()) {
    return false;
  }

  cam = BrownPinholeCamera(f, ppx, ppy, R, t, k1, k2, k3);
  return true;
}

/// Save a BrownPinholeCamera to a ascii file
static bool save(
  const std::string & scameraFile,
  const BrownPinholeCamera & cam)
{
  if (stlplus::extension_part(scameraFile) != "txt")
    return false;

  std::ofstream file(scameraFile.c_str(), std::ios::out);
  // Save intrinsic data:
  file << cam._f << " "
    << cam._ppx << " "
    << cam._ppy << " "
    << cam._k1 << " "
    << cam._k2 << " "
    << cam._k3 << "\n";
  // Save extrinsic data
  const Mat3 & R = cam._R;
  file << R(0,0) << " " << R(0,1) << " " << R(0,2) << "\n"
    << R(1,0) << " " << R(1,1) << " " << R(1,2) << "\n"
    << R(2,0) << " " << R(2,1) << " " << R(2,2) << "\n";
  file << cam._t(0) << " " << cam._t(1) << " " << cam._t(2) << "\n";
  bool bOk = (!file.fail());
  file.close();
  return bOk;
}

} // namespace openMVG

#endif // OPENMVG_CAMERA_IO_H
