
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_IO_H
#define OPENMVG_CAMERA_IO_H

#include "openMVG/cameras/PinholeCamera.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>

namespace openMVG
{
namespace cameras
{


/**
* @brief Save a Pinhole camera to a file as a Projection matrix (12 doubles as binary values)
* @param scameraFile path to the file in which matrix will be written
* @param cam Camera to write
* @retval true If write is correct
* @retval false If there was an error during export
*/
static bool save(
  const std::string & scameraFile,
  const PinholeCamera & cam )
{
  if ( stlplus::extension_part( scameraFile ) != "bin" )
  {
    return false;
  }

  const Mat34 & PMat = cam._P;
  std::ofstream file( scameraFile.c_str(), std::ios::out | std::ios::binary );
  file.write( ( const char* )PMat.data(), ( std::streamsize )( 3 * 4 )*sizeof( double ) );

  const bool bOk = ( !file.fail() );
  file.close();
  return bOk;
}


/**
* @brief Load a Pinhole camera from a file (read 12 doubles saved in binary values)
* @param scameraFile path where camera is stored
* @param[out] cam output camera after read
* @retval true If loading is correct
* @retval false If there was an error during loading
*/
static bool load(
  const std::string & scameraFile,
  PinholeCamera & cam )
{
  std::vector<double> val;

  if ( stlplus::extension_part( scameraFile ) == "bin" )
  {
    std::ifstream in( scameraFile.c_str(), std::ios::in | std::ios::binary );
    if ( !in.is_open() )
    {
      std::cerr << "Error: failed to open file '" << scameraFile
                << "' for reading" << std::endl;
      return false;
    }
    val.resize(12);
    in.read(reinterpret_cast<char*>(&val[0]),(std::streamsize)12*sizeof(double));
    if (in.fail())
    {
      val.clear();
    }
  }
  else
  {
    return false;
  }

  if ( val.size() == 3 * 4 ) //P Matrix
  {
    Mat34 P;
    P << val[0], val[3], val[6], val[9],
    val[1], val[4], val[7], val[10],
    val[2], val[5], val[8], val[11];

    Mat3 R, K;
    Vec3 t;
    KRt_From_P( P, &K, &R, &t );
    cam = PinholeCamera( K, R, t );
    return true;
  }
  return false;
}

} // namespace cameras
} // namespace openMVG

#endif // OPENMVG_CAMERA_IO_H
