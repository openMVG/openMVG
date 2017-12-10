// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <cmath>
#include <iterator>
#include <iomanip>
#include <fstream>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;


bool exportToPMVSFormat(
  const SfM_Data & sfm_data,
  const std::string & sOutDirectory,  //Output PMVS files directory
  const int downsampling_factor,
  const int CPU_core_count,
  const bool b_VisData = true
  )
{
  bool bOk = true;
  if (!stlplus::is_folder(sOutDirectory))
  {
    stlplus::folder_create(sOutDirectory);
    bOk = stlplus::is_folder(sOutDirectory);
  }

  // Create basis directory structure
  stlplus::folder_create( stlplus::folder_append_separator(sOutDirectory) + "models");
  stlplus::folder_create( stlplus::folder_append_separator(sOutDirectory) + "txt");
  stlplus::folder_create( stlplus::folder_append_separator(sOutDirectory) + "visualize");

  if (bOk &&
      stlplus::is_folder(stlplus::folder_append_separator(sOutDirectory) + "models") &&
      stlplus::is_folder( stlplus::folder_append_separator(sOutDirectory) + "txt") &&
      stlplus::is_folder( stlplus::folder_append_separator(sOutDirectory) + "visualize")
      )
  {
    bOk = true;
  }
  else  {
    std::cerr << "Cannot access to one of the desired output directory" << std::endl;
  }

  if (bOk)
  {
    C_Progress_display my_progress_bar( sfm_data.GetViews().size()*2 );

    // Since PMVS requires contiguous camera index, and that some views can have some missing poses,
    // we reindex the poses to ensure a contiguous pose list.
    Hash_Map<IndexT, IndexT> map_viewIdToContiguous;

    // Export valid views as Projective Cameras:
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
        iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
    {
      const View * view = iter->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      // View Id re-indexing
      map_viewIdToContiguous.insert(std::make_pair(view->id_view, map_viewIdToContiguous.size()));

      // We have a valid view with a corresponding camera & pose
      const Mat34 P = iterIntrinsic->second->get_projective_equivalent(pose);
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << map_viewIdToContiguous[view->id_view];
      std::ofstream file(
        stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory) + "txt",
        os.str() ,"txt").c_str());
      file << "CONTOUR" << os.widen('\n')
        << P.row(0) <<"\n"<< P.row(1) <<"\n"<< P.row(2) << os.widen('\n');
      file.close();
    }

    // Export (calibrated) views as undistorted images
    Image<RGBColor> image, image_ud;
    const Views & views = sfm_data.GetViews();
    #pragma omp parallel for private(image, image_ud)
    for (int i = 0; i < static_cast<int>(views.size()); ++i)
    {
      ++my_progress_bar;
      auto view_it = views.begin();
      std::advance(view_it, i);
      const View * view = view_it->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      // We have a valid view with a corresponding camera & pose
      const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << map_viewIdToContiguous[view->id_view];
      const std::string dstImage = stlplus::create_filespec(
        stlplus::folder_append_separator(sOutDirectory) + "visualize", os.str(),"jpg");

      const IntrinsicBase * cam = iterIntrinsic->second.get();
      if (cam->have_disto())
      {
        // undistort the image and save it
        ReadImage( srcImage.c_str(), &image);
        UndistortImage(image, cam, image_ud, BLACK);
        WriteImage(dstImage.c_str(), image_ud);
      }
      else // (no distortion)
      {
        // copy the image if extension match
        if (stlplus::extension_part(srcImage) == "JPG" ||
          stlplus::extension_part(srcImage) == "jpg")
        {
          stlplus::file_copy(srcImage, dstImage);
        }
        else
        {
          ReadImage( srcImage.c_str(), &image);
          WriteImage( dstImage.c_str(), image);
        }
      }
    }

    //pmvs_options.txt
    std::ostringstream os;
    os << "level " << downsampling_factor << os.widen('\n')
     << "csize 2" << os.widen('\n')
     << "threshold 0.7" << os.widen('\n')
     << "wsize 7" << os.widen('\n')
     << "minImageNum 3" << os.widen('\n')
     << "CPU " << CPU_core_count << os.widen('\n')
     << "setEdge 0" << os.widen('\n')
     << "useBound 0" << os.widen('\n')
     << "useVisData " << (int) b_VisData << os.widen('\n')
     << "sequence -1" << os.widen('\n')
     << "maxAngle 10" << os.widen('\n')
     << "quad 2.0" << os.widen('\n')
     << "timages -1 0 " << map_viewIdToContiguous.size() << os.widen('\n')
     << "oimages 0" << os.widen('\n');

    if (b_VisData)
    {
      std::map<IndexT, std::set<IndexT>> view_shared;
      // From the structure observations, list the putatives pairs (symmetric)
      for (Landmarks::const_iterator itL = sfm_data.GetLandmarks().begin();
        itL != sfm_data.GetLandmarks().end(); ++itL)
      {
        const Landmark & landmark = itL->second;
        const Observations & obs = landmark.obs;
        for (Observations::const_iterator itOb = obs.begin();
          itOb != obs.end(); ++itOb)
        {
          const IndexT viewId = itOb->first;
          Observations::const_iterator itOb2 = itOb;
          ++itOb2;
          for (; itOb2 != obs.end(); ++itOb2)
          {
            const IndexT viewId2 = itOb2->first;
            view_shared[map_viewIdToContiguous[viewId]].insert(map_viewIdToContiguous[viewId2]);
            view_shared[map_viewIdToContiguous[viewId2]].insert(map_viewIdToContiguous[viewId]);
          }
        }
      }
      // Export the vis.dat file
      std::ostringstream osVisData;
      osVisData
        << "VISDATA" << os.widen('\n')
        << view_shared.size() << os.widen('\n'); // #images
      // Export view shared visibility
      for (std::map<IndexT, std::set<IndexT>>::const_iterator it = view_shared.begin();
        it != view_shared.end(); ++it)
      {
        const std::set<IndexT> & setView = it->second;
        osVisData << it->first << ' ' << setView.size();
        for (std::set<IndexT>::const_iterator itV = setView.begin();
          itV != setView.end(); ++itV)
        {
          osVisData << ' ' << *itV;
        }
        osVisData << os.widen('\n');
      }
      std::ofstream file(stlplus::create_filespec(sOutDirectory, "vis", "dat").c_str());
      file << osVisData.str();
      file.close();
    }

    std::ofstream file(stlplus::create_filespec(sOutDirectory, "pmvs_options", "txt").c_str());
    file << os.str();
    file.close();
  }
  return bOk;
}

bool exportToBundlerFormat(
  const SfM_Data & sfm_data,
  const std::string & sOutFile, //Output Bundle.rd.out file
  const std::string & sOutListFile)  //Output Bundler list.txt file
{
  std::ofstream os(sOutFile.c_str());
  std::ofstream osList(sOutListFile.c_str());
  if (! os.is_open() || ! osList.is_open())
  {
    return false;
  }
  else
  {
    // Since PMVS requires contiguous camera index, and that some views can have some missing poses,
    // we reindex the poses to ensure a contiguous pose list.
    Hash_Map<IndexT, IndexT> map_viewIdToContiguous;

    // Count the number of valid cameras and re-index the viewIds
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
        iter != sfm_data.GetViews().end(); ++iter)
    {
      const View * view = iter->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      // View Id re-indexing
      map_viewIdToContiguous.insert(std::make_pair(view->id_view, map_viewIdToContiguous.size()));
    }

    // Fill the "Bundle file"
    os << "# Bundle file v0.3" << os.widen('\n')
      << map_viewIdToContiguous.size()  << " " << sfm_data.GetLandmarks().size() << os.widen('\n');

    // Export camera properties & image filenames
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
        iter != sfm_data.GetViews().end(); ++iter)
    {
      const View * view = iter->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      // Must export focal, k1, k2, R, t

      Mat3 D;
      D.fill(0.0);
      D .diagonal() = Vec3(1., -1., -1.); // mapping between our pinhole and Bundler convention
      const double k1 = 0.0, k2 = 0.0; // distortion already removed

      if (isPinhole(iterIntrinsic->second->getType()))
      {
        const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic->second.get());
        const double focal = cam->focal();
        const Mat3 R = D * pose.rotation();
        const Vec3 t = D * pose.translation();

        os << focal << " " << k1 << " " << k2 << os.widen('\n') //f k1 k2
          << R(0,0) << " " << R(0, 1) << " " << R(0, 2) << os.widen('\n')  //R.row(0)
          << R(1,0) << " " << R(1, 1) << " " << R(1, 2) << os.widen('\n')  //R.row(1)
          << R(2,0) << " " << R(2, 1) << " " << R(2, 2) << os.widen('\n')  //R.row(2)
          << t(0)   << " " << t(1)    << " " << t(2)    << os.widen('\n'); //t

        osList << stlplus::basename_part(view->s_Img_path) + "." + stlplus::extension_part(view->s_Img_path)
          << " 0 " << focal << os.widen('\n');
      }
      else
      {
        std::cerr << "Unsupported camera model for Bundler export." << std::endl;
        return false;
      }
    }
    // Export structure and visibility
    for (Landmarks::const_iterator iter = sfm_data.GetLandmarks().begin();
      iter != sfm_data.GetLandmarks().end(); ++iter)
    {
      const Landmark & landmark = iter->second;
      const Observations & obs = landmark.obs;
      const Vec3 & X = landmark.X;
      // X, color, obsCount
      os << X[0] << " " << X[1] << " " << X[2] << os.widen('\n')
         <<  "255 255 255" << os.widen('\n')
         << obs.size() << " ";
      for (Observations::const_iterator iterObs = obs.begin();
        iterObs != obs.end(); ++iterObs)
      {
        const Observation & ob = iterObs->second;
        // ViewId, FeatId, x, y
        os << map_viewIdToContiguous[iterObs->first] << " " << ob.id_feat << " " << ob.x(0) << " " << ob.x(1) << " ";
      }
      os << os.widen('\n');
    }
    os.close();
    osList.close();
  }
  return true;
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  int resolution = 1;
  int CPU = 8;
  bool bVisData = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', resolution, "resolution") );
  cmd.add( make_option('c', CPU, "CPU") );
  cmd.add( make_option('v', bVisData, "useVisData") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir path]\n"
      << "[-r|--resolution] divide image coefficient\n"
      << "[-c|--nb core]\n"
      << "[-v|--useVisData] use visibility information."
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  {
    exportToPMVSFormat(sfm_data,
      stlplus::folder_append_separator(sOutDir) + "PMVS",
      resolution,
      CPU,
      bVisData);

    exportToBundlerFormat(sfm_data,
      stlplus::folder_append_separator(sOutDir) +
      stlplus::folder_append_separator("PMVS") + "bundle.rd.out",
      stlplus::folder_append_separator(sOutDir) +
      stlplus::folder_append_separator("PMVS") + "list.txt"
      );

    return EXIT_SUCCESS;
  }

  // Exit program
  return EXIT_FAILURE;
}
