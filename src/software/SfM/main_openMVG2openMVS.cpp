// Copyright (c) 2016 
// cDc <cdc.seacave@gmail.com>
// Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image.hpp"

#define _USE_EIGEN
#include "InterfaceMVS.h"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>
#include <iomanip>

bool exportToOpenMVS(
  const SfM_Data & sfm_data,
  const std::string & sOutFile,
  const std::string & sOutDir
  )
{
  // Create undistorted images directory structure
  if (!stlplus::is_folder(sOutDir))
  {
    stlplus::folder_create(sOutDir);
    if (!stlplus::is_folder(sOutDir))
    {
      std::cerr << "Cannot access to one of the desired output directory" << std::endl;
      return false;
    }
  }

  // Export data :
  MVS::Interface scene;
  size_t nPoses(0);
  const size_t nCameras(sfm_data.GetIntrinsics().size());
  const uint32_t nViews((uint32_t)sfm_data.GetViews().size());

  C_Progress_display my_progress_bar(nViews);

  // OpenMVG can have not contiguous index, use a map to create the required OpenMVS contiguous ID index
  std::map<openMVG::IndexT, uint32_t> map_intrinsic, map_view;

  // define a platform with all the intrinsic group
  for (const auto& intrinsic: sfm_data.GetIntrinsics())
  {
    if (isPinhole(intrinsic.second.get()->getType()))
    {
      const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(intrinsic.second.get());
      if (map_intrinsic.count(intrinsic.first) == 0)
        map_intrinsic.insert(std::make_pair(intrinsic.first, scene.platforms.size()));
      MVS::Interface::Platform platform;
      // add the camera
      MVS::Interface::Platform::Camera camera;
      camera.K = cam->K();
      // sub-pose
      camera.R = Mat3::Identity();
      camera.C = Vec3::Zero();
      platform.cameras.push_back(camera);
      scene.platforms.push_back(platform);
    }
  }

  // define images & poses
  scene.images.reserve(nViews);
  for (const auto& view : sfm_data.GetViews())
  {
    map_view[view.first] = scene.images.size();
    MVS::Interface::Image image;
    const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view.second->s_Img_path);
    image.name = stlplus::create_filespec(sOutDir, view.second->s_Img_path);
    image.platformID = map_intrinsic.at(view.second->id_intrinsic);
    MVS::Interface::Platform& platform = scene.platforms[image.platformID];
    image.cameraID = 0;
    if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get()) && stlplus::is_file(srcImage))
    {
      MVS::Interface::Platform::Pose pose;
      image.poseID = platform.poses.size();
      const openMVG::geometry::Pose3 poseMVG(sfm_data.GetPoseOrDie(view.second.get()));
      pose.R = poseMVG.rotation();
      pose.C = poseMVG.center();
      // export undistorted images
      const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view.second->id_intrinsic).get();
      if (cam->have_disto())
      {
        // undistort image and save it
        Image<openMVG::image::RGBColor> imageRGB, imageRGB_ud;
        ReadImage(srcImage.c_str(), &imageRGB);
        UndistortImage(imageRGB, cam, imageRGB_ud, BLACK);
        WriteImage(image.name.c_str(), imageRGB_ud);
      }
      else
      {
        // just copy image
        stlplus::file_copy(srcImage, image.name);
      }
      platform.poses.push_back(pose);
      ++nPoses;
    }
    else
    {
      // image have not valid pose, so set an undefined pose
      image.poseID = NO_ID;
      // just copy the image
      stlplus::file_copy(srcImage, image.name);
    }
    scene.images.push_back(image);
    ++my_progress_bar;
  }

  // define structure
  scene.vertices.reserve(sfm_data.GetLandmarks().size());
  for (const auto& vertex: sfm_data.GetLandmarks())
  {
    const Landmark & landmark = vertex.second;
    MVS::Interface::Vertex vert;
    MVS::Interface::Vertex::ViewArr& views = vert.views;
    for (const auto& observation: landmark.obs)
    {
      const auto it(map_view.find(observation.first));
      if (it != map_view.end()) {
        MVS::Interface::Vertex::View view;
        view.imageID = it->second;
        view.confidence = 0;
        views.push_back(view);
      }
    }
    if (views.size() < 2)
      continue;
    std::sort(
      views.begin(), views.end(),
      [] (const MVS::Interface::Vertex::View& view0, const MVS::Interface::Vertex::View& view1)
      {
        return view0.imageID < view1.imageID;
      }
    );
    vert.X = landmark.X.cast<float>();
    scene.vertices.push_back(vert);
  }

  // normalize camera intrinsics
  for (size_t p=0; p<scene.platforms.size(); ++p)
  {
    MVS::Interface::Platform& platform = scene.platforms[p];
    for (size_t c=0; c<platform.cameras.size(); ++c) {
      MVS::Interface::Platform::Camera& camera = platform.cameras[c];
      // find one image using this camera
      MVS::Interface::Image* pImage(NULL);
      for (MVS::Interface::Image& image: scene.images)
      {
	      if (image.platformID == p && image.cameraID == c && image.poseID != NO_ID)
	      {
		      pImage = &image;
		      break;
	      }
      }
      if (pImage == NULL)
      {
	      std::cerr << "error: no image using camera " << c << " of platform " << p << std::endl;
	      continue;
      }
      // read image meta-data
      ImageHeader imageHeader;
      ReadImageHeader(pImage->name.c_str(), &imageHeader);
      const double fScale(1.0/std::max(imageHeader.width, imageHeader.height));
      camera.K(0, 0) *= fScale;
      camera.K(1, 1) *= fScale;
      camera.K(0, 2) *= fScale;
      camera.K(1, 2) *= fScale;
    }
  }

  // write OpenMVS data
  if (!ARCHIVE::SerializeSave(scene, sOutFile))
    return false;

  std::cout
    << "Scene saved to OpenMVS interface format:\n"
    << "\t" << scene.images.size() << " images (" << nPoses << " calibrated)\n"
    << "\t" << scene.vertices.size() << " Landmarks\n";
  return true;
}

int main(int argc, char *argv[])
{
  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutFile = "scene.mvs";
  std::string sOutDir = "undistorted_images";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutFile, "outfile") );
  cmd.add( make_option('d', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outfile] OpenMVS scene file\n"
      << "[-d|--outdir] undistorted images path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (exportToOpenMVS(sfm_data, sOutFile, sOutDir))
    return( EXIT_SUCCESS );
  return( EXIT_FAILURE );
}
