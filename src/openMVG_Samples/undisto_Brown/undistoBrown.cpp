
// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/numeric/numeric.h"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <string>
#include <iostream>

using namespace openMVG;
using namespace std;

// A simple container and undistort function for the Brown's distortion model [1]
// Variables:
// (x,y): 2D point in the image (pixel)
// (u,v): the undistorted 2D point (pixel)
// radial_distortion (k1, k2, k3, ...): vector containing the radial distortion
// (cx,cy): camera principal point
// tangential factors are not considered here
//
// Equation:
// u = x + (x - cx) * (k1 * r^2 + k2 * r^4 +...)
// v = y + (y - cy) * (k1 * r^2 + k2 * r^4 +...)
//
// [1] Decentering distortion of lenses.
//      Brown, Duane C
//      Photometric Engineering 1966
struct BrownDistoModel
{
  Vec2 m_disto_center; // distortion center
  Vec m_radial_distortion; // radial distortion factor
  double m_f; // focal

  inline void ComputeUndistortedCoordinates( double xu, double yu, double &xd, double& yd) const
  {
    const Vec2 point (xu, yu);
    const Vec2 principal_point (m_disto_center);
    const Vec2 point_centered = point - principal_point;

    const double u = point_centered.x() / m_f;
    const double v = point_centered.y() / m_f;
    const double radius_squared = u * u + v * v;

    double coef_radial = 0.0;
    for (int i = m_radial_distortion.size() - 1; i >= 0; --i) {
      coef_radial = (coef_radial + m_radial_distortion[i]) * radius_squared;
    }

    const Vec2 undistorted_point = point + point_centered * coef_radial;
    xd = undistorted_point(0);
    yd = undistorted_point(1);
  }
};

/// Undistort an image according a given Distortion model
template <typename Image>
Image undistortImage(
  const Image& I,
  const BrownDistoModel& d,
  RGBColor fillcolor = BLACK,
  bool bcenteringPPpoint = true)
{
  const int w = I.Width();
  const int h = I.Height();
  Vec2 offset(0,0);
  if (bcenteringPPpoint)
    offset = Vec2(w * .5, h * .5) - d.m_disto_center;

  Image J ( w,h );
  double xu, yu, xd,yd;
  for ( int j=0; j<h; j++ ) {
    for ( int i=0; i<w; i++ ) {
      xu = double (i);
      yu = double (j);
      d.ComputeUndistortedCoordinates(xu, yu, xd, yd);
      xd -= offset(0);
      yd -= offset(1);
      if ( !J.Contains(yd, xd) )
        J ( j,i ) = fillcolor;
      else
        J ( j,i ) = SampleLinear(I, yd, xd);
    }
  }
  return J;
}

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sPath;
  std::string sOutPath;
  // Temp storage for the Brown's distortion model
  Vec2 c; // distortion center
  Vec3 k; // distortion factor
  double f; // Focal
  std::string suffix = "JPG";

  cmd.add( make_option('i', sPath, "imadir") );
  cmd.add( make_option('o', sOutPath, "outdir") );
  cmd.add( make_option('a', c(0), "cx") );
  cmd.add( make_option('b', c(1), "cy") );
  cmd.add( make_option('c', k(0), "k1") );
  cmd.add( make_option('d', k(1), "k2") );
  cmd.add( make_option('e', k(2), "k3") );
  cmd.add( make_option('f', f, "focal") );
  cmd.add( make_option('s', suffix, "suffix") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--imadir - Input path]"
      << "[-o|--outdir - path for the undistorted JPG files]"
      << "[-f|--focal - focal length]"
      << "[-s|--suffix - Suffix of the input files. (default: JPG)]"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutPath == sPath)
  {
    std::cerr << "Input and Ouput path are set to the same value" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutPath))
    stlplus::folder_create(sOutPath);

  BrownDistoModel distoModel;
  distoModel.m_disto_center = Vec2(c(0), c(1));
  distoModel.m_radial_distortion = k;
  distoModel.m_f = f;

  std::cout << "Used Brown's distortion model values: \n"
    << "  Distortion center: " << distoModel.m_disto_center.transpose() << "\n"
    << "  Distortion coefficients (K1,K2,K3): "
    << distoModel.m_radial_distortion.transpose() << "\n"
    << "  Distortion focal: " << distoModel.m_f << std::endl;

  const std::vector<std::string> vec_fileNames =
    stlplus::folder_wildcard(sPath, "*."+suffix, false, true);
  std::cout << "\nLocated " << vec_fileNames.size() << " files in " << sPath
    << " with suffix " << suffix;

  Image<unsigned char > imageGreyIn, imageGreyU;
  Image<RGBColor> imageRGBIn, imageRGBU;
  Image<RGBAColor> imageRGBAIn, imageRGBAU;
  
  C_Progress_display my_progress_bar( vec_fileNames.size() );
  for (size_t j = 0; j < vec_fileNames.size(); ++j, ++my_progress_bar)
  {
    //read the depth
    int w,h,depth;
    vector<unsigned char> tmp_vec; 
    const string sOutFileName =
      stlplus::create_filespec(sOutPath, stlplus::basename_part(vec_fileNames[j]), "JPG");
    const string sInFileName = stlplus::create_filespec(sPath, stlplus::basename_part(vec_fileNames[j]));
    const int res = ReadImage(sInFileName.c_str(), &tmp_vec, &w, &h, &depth);
    if (res == 1)
    {
      switch(depth)
      {
        case 1: //Greyscale
          {
            imageGreyIn = Eigen::Map<Image<unsigned char>::Base>(&tmp_vec[0], h, w);
            imageGreyU = undistortImage ( imageGreyIn, distoModel);
            WriteImage(sOutFileName.c_str(), imageGreyU);
            break;
          }
        case 3: //RGB
          {
            imageRGBIn = Eigen::Map<Image<RGBColor>::Base>((RGBColor*) &tmp_vec[0], h, w);
            imageRGBU = undistortImage ( imageRGBIn, distoModel);
            WriteImage(sOutFileName.c_str(), imageRGBU);
            break;
          }
        case 4: //RGBA
          {
            imageRGBAIn = Eigen::Map<Image<RGBAColor>::Base>((RGBAColor*) &tmp_vec[0], h, w);
            imageRGBAU = undistortImage ( imageRGBAIn, distoModel);
            WriteImage(sOutFileName.c_str(), imageRGBAU);
            break;
          }
      }

    }//end if res==1
    else
    {
      std::cerr << "\nThe image contains " << depth << "layers. This depth is not supported!\n";
    }
  } //end loop for each file
  return EXIT_SUCCESS;
}

