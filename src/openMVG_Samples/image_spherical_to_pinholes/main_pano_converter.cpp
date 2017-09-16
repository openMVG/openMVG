// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/sample.hpp"
#include "./panorama_helper.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>

using namespace std;
using namespace svg;

// Convert spherical panorama to rectilinear images
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string
    s_directory_in,
    s_directory_out;
  int
    image_resolution = 1200,
    nb_split = 5;

  // required
  cmd.add( make_option('i', s_directory_in, "input_dir") );
  cmd.add( make_option('o', s_directory_out, "output_dir") );
  // Optional
  cmd.add( make_option('r', image_resolution, "image_resolution") );
  cmd.add( make_option('n', nb_split, "nb_split") );
  cmd.add( make_switch('D', "demo_mode") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_dir] the path where the spherical panoramic images are saved \n"
    << "[-o|--output_dir] the path where output rectilinear image will be saved \n"
    << " OPTIONAL:\n"
    << "[-r|--image_resolution] the rectilinear image size (default:" << image_resolution << ") \n"
    << "[-n|--nb_split] the number of rectilinear image along the X axis (default:" << nb_split << ") \n"
    << "[-D|--demo_mode] switch parameter, export a SVG file that simulate asked rectilinear\n"
    << "  frustum configuration on the spherical image.\n"
    << std::endl;
  }

  // Input parameter checking

  if (image_resolution < 0)
  {
    std::cerr << "image_resolution must be larger than 0" << std::endl;
    return EXIT_FAILURE;
  }
  if (nb_split < 0)
  {
    std::cerr << "nb_split must be larger than 0" << std::endl;
    return EXIT_FAILURE;
  }
  if (s_directory_in.empty() || s_directory_out.empty())
  {
    std::cerr << "input_dir and output_dir option must not be empty" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::is_folder(s_directory_out))
  {
    if (!stlplus::folder_create(s_directory_out))
    {
      std::cerr << "Cannot create the output_dir directory" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // List images from the input directory
  const std::vector<std::string> vec_filenames
    = stlplus::folder_wildcard(s_directory_in, "*.jpg", false, true);

  if (vec_filenames.empty())
  {
    std::cerr << "Did not find any jpg image in the provided input_dir" << std::endl;
    return EXIT_FAILURE;
  }

  using namespace openMVG;

  //-----------------
  //-- Create N rectilinear cameras
  //     (according the number of asked split along the X axis)
  //-- For each camera
  //   -- Forward mapping
  //   -- Save resulting images to disk
  //-----------------

  using CGeomFunctor = CsphericalMapping;

  //-- Generate N cameras along the X axis
  std::vector< openMVG::PinholeCamera_R > vec_cam;

  const double twoPi = M_PI * 2.0;
  const double alpha = twoPi / static_cast<double>(nb_split);

  const int wIma = image_resolution, hIma = image_resolution;
  const double focal = openMVG::focalFromPinholeHeight(hIma, openMVG::D2R(60));
  double angle = 0.0;
  for (int i = 0; i < nb_split; ++i, angle += alpha)
  {
    vec_cam.emplace_back(focal, wIma, hIma, RotationAroundY(angle));
  }

  if (cmd.used('D')) // Demo mode:
  {
    const int wPano = 4096, hPano = wPano / 2;

    svgDrawer svgStream(wPano, hPano);
    svgStream.drawLine(0,0,wPano,hPano, svgStyle());
    svgStream.drawLine(wPano,0, 0, hPano, svgStyle());

    //--> For each cam, reproject the image borders onto the panoramic image

    for (const openMVG::PinholeCamera_R & cam_it : vec_cam)
    {
      //draw the shot border with the givenStep:
      const int step = 10;
      Vec3 ray;

      // Vertical rectilinear image border:
      for (double j = 0; j <= image_resolution; j += image_resolution/(double)step)
      {
        Vec2 pt(0.,j);
        ray = cam_it.getRay(pt(0), pt(1));
        Vec2 x = CGeomFunctor::Get2DPoint( ray, wPano, hPano);
        svgStream.drawCircle(x(0), x(1), 4, svgStyle().fill("green"));

        pt[0] = image_resolution;
        ray = cam_it.getRay(pt(0), pt(1));
        x = CGeomFunctor::Get2DPoint( ray, wPano, hPano);
        svgStream.drawCircle(x(0), x(1), 4, svgStyle().fill("green"));
      }
      // Horizontal rectilinear image border:
      for (double j = 0; j <= image_resolution; j += image_resolution/(double)step)
      {
        Vec2 pt(j,0.);
        ray = cam_it.getRay(pt(0), pt(1));
        Vec2 x = CGeomFunctor::Get2DPoint( ray, wPano, hPano);
        svgStream.drawCircle(x(0), x(1), 4, svgStyle().fill("yellow"));

        pt[1] = image_resolution;
        ray = cam_it.getRay(pt(0), pt(1));
        x = CGeomFunctor::Get2DPoint( ray, wPano, hPano);
        svgStream.drawCircle(x(0), x(1), 4, svgStyle().fill("yellow"));
      }
    }

    std::ofstream svgFile( stlplus::create_filespec(s_directory_out, "test.svg" ));
    svgFile << svgStream.closeSvgFile().str();

    return EXIT_SUCCESS;
  }


  //-- For each input image extract multiple pinhole images
  for (const std::string & filename_it : vec_filenames)
  {
    image::Image<image::RGBColor> imageSource;
    if (!ReadImage(stlplus::create_filespec(s_directory_in,filename_it).c_str(), &imageSource))
    {
      std::cerr << "Cannot read the image" << std::endl;
      continue;
    }

    const int
      wPano = imageSource.Width(),
      hPano = imageSource.Height();

    const image::Sampler2d<image::SamplerLinear> sampler;
    image::Image<image::RGBColor> imaOut(wIma, hIma, image::BLACK);

    size_t index = 0;
    for (const PinholeCamera_R & cam_it : vec_cam)
    {
      imaOut.fill(image::BLACK);

      // Backward mapping:
      // - Find for each pixels of the pinhole image where it comes from the panoramic image
      for (int j = 0; j < hIma; ++j)
      {
        for (int i = 0; i < wIma; ++i)
        {
          const Vec3 ray = cam_it.getRay(i, j);
          const Vec2 x = CGeomFunctor::Get2DPoint(ray, wPano, hPano);
          imaOut(j,i) = sampler(imageSource, x(1), x(0));
        }
      }
      //-- save image
      const std::string basename = stlplus::basename_part(filename_it);

      std::cout << basename << " cam index: " << index << std::endl;

      std::ostringstream os;
      os << s_directory_out << "/" << basename << "_" << index << ".jpg";
      WriteImage(os.str().c_str(), imaOut);

      ++index;
    }
  }

  std::ofstream fileFocalOut(stlplus::create_filespec(s_directory_out, "focal.txt"));
  fileFocalOut << focal;
  fileFocalOut.close();

  return EXIT_SUCCESS;
}
