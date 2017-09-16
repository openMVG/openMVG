
// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/image/image_io.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::image;
using namespace std;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sPath;
  std::string sOutPath;
  // Temp storage for the Brown's distortion model
  Vec2 c; // distortion center
  Vec3 k; // distortion factors
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
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--imadir - Input path]\n"
      << "[-o|--outdir - path for the undistorted JPG files]\n"
      << "[-f|--focal - focal length]\n"
      << "[-s|--suffix - Suffix of the input files. (default: JPG)]\n"
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

  std::cout << "Used Brown's distortion model values: \n"
    << "  Distortion center: " << c.transpose() << "\n"
    << "  Distortion coefficients (K1,K2,K3): "
    << k.transpose() << "\n"
    << "  Distortion focal: " << f << std::endl;

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
      stlplus::create_filespec(sOutPath, stlplus::basename_part(vec_fileNames[j]), "png");
    const string sInFileName = stlplus::create_filespec(sPath, stlplus::filename_part(vec_fileNames[j]));
    const int res = ReadImage(sInFileName.c_str(), &tmp_vec, &w, &h, &depth);

    const Pinhole_Intrinsic_Radial_K3 cam(w, h, f, c(0), c(1), k(0), k(1), k(2));

    if (res == 1)
    {
      switch (depth)
      {
        case 1: //Greyscale
          {
            imageGreyIn = Eigen::Map<Image<unsigned char>::Base>(&tmp_vec[0], h, w);
            UndistortImage(imageGreyIn, &cam, imageGreyU);
            WriteImage(sOutFileName.c_str(), imageGreyU);
            break;
          }
        case 3: //RGB
          {
            imageRGBIn = Eigen::Map<Image<RGBColor>::Base>((RGBColor*) &tmp_vec[0], h, w);
            UndistortImage(imageRGBIn, &cam, imageRGBU);
            WriteImage(sOutFileName.c_str(), imageRGBU);
            break;
          }
        case 4: //RGBA
          {
            imageRGBAIn = Eigen::Map<Image<RGBAColor>::Base>((RGBAColor*) &tmp_vec[0], h, w);
            UndistortImage(imageRGBAIn, &cam, imageRGBAU);
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
