// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/image/image_io.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"

#include <cstdlib>
#include <cmath>
#include <iterator>
#include <iomanip>
#include <fstream>

bool exportToCMPMVSFormat(
  const SfM_Data & sfm_data,
  const std::string & sOutDirectory // Output CMPMVS files directory
  )
{
  bool bOk = true;
  // Create basis directory structure
  if (!stlplus::is_folder(sOutDirectory))
  {
    stlplus::folder_create(sOutDirectory);
    bOk = stlplus::is_folder(sOutDirectory);
  }

  if (!bOk)
  {
    std::cerr << "Cannot access to one of the desired output directory" << std::endl;
    return false;
  }
  else
  {
    // Export data:

    C_Progress_display my_progress_bar( sfm_data.GetViews().size()*2 );

    // Since CMPMVS requires contiguous camera index, and that some views can have some missing poses,
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
      os << std::setw(5) << std::setfill('0') << map_viewIdToContiguous[view->id_view] << "_P";
      std::ofstream file(
        stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory),
        os.str() ,"txt").c_str());
      file << "CONTOUR" << os.widen('\n')
        << P.row(0) <<"\n"<< P.row(1) <<"\n"<< P.row(2) << os.widen('\n');
      file.close();
    }

    // Export (calibrated) views as undistorted images
    std::pair<unsigned int, unsigned int> w_h_image_size;
    Image<RGBColor> image, image_ud;
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
        iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
    {
      const View * view = iter->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);

      // We have a valid view with a corresponding camera & pose
      const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
      std::ostringstream os;
      os << std::setw(5) << std::setfill('0') << map_viewIdToContiguous[view->id_view];
      std::string dstImage = stlplus::create_filespec(
        stlplus::folder_append_separator(sOutDirectory), os.str(),"jpg");

      const IntrinsicBase * cam = iterIntrinsic->second.get();
      if (map_viewIdToContiguous[view->id_view] == 0)
        w_h_image_size = std::make_pair(cam->w(), cam->h());
      else
      {
        // check that there is no image sizing change (CMPMVS support only images of the same size)
        if (cam->w() != w_h_image_size.first ||
            cam->h() != w_h_image_size.second)
        {
          std::cerr << "CMPMVS support only image having the same image size";
          return false;
        }
      }
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

    // Write the mvs_firstRun script
    std::ostringstream os;
    os << "[global]" << os.widen('\n')
      << "dirName=\"" << stlplus::folder_append_separator(sOutDirectory) <<"\"" << os.widen('\n')
      << "prefix=\"\"" << os.widen('\n')
      << "imgExt=\"jpg\"" << os.widen('\n')
      << "ncams=" << map_viewIdToContiguous.size() << os.widen('\n')
      << "width=" << w_h_image_size.first << os.widen('\n')
      << "height=" << w_h_image_size.second << os.widen('\n')
      << "scale=2" << os.widen('\n')
      << "workDirName=\"_tmp_fast\"" << os.widen('\n')
      << "doPrepareData=TRUE" << os.widen('\n')
      << "doPrematchSifts=TRUE" << os.widen('\n')
      << "doPlaneSweepingSGM=TRUE"  << os.widen('\n')
      << "doFuse=TRUE" << os.widen('\n')
      << "nTimesSimplify=10" << os.widen('\n')
      << os.widen('\n')
      << "[prematching]" << os.widen('\n')
      << "minAngle=3.0" << os.widen('\n')
      << os.widen('\n')
      << "[grow]" << os.widen('\n')
      << "minNumOfConsistentCams=6" << os.widen('\n')
      << os.widen('\n')
      << "[filter]" << os.widen('\n')
      << "minNumOfConsistentCams=2" << os.widen('\n')
      << os.widen('\n')
      << "#do not erase empy lines after this comment otherwise it will crash ... bug" << os.widen('\n')
      << os.widen('\n')
      << os.widen('\n');

    std::ofstream file(
      stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory),
      "01_mvs_firstRun" ,"ini").c_str());
    file << os.str();
    file.close();

    // limitedScale
    os.str("");
    os << "[global]" << os.widen('\n')
      << "dirName=\"" << stlplus::folder_append_separator(sOutDirectory) <<"\"" << os.widen('\n')
      << "prefix=\"\"" << os.widen('\n')
      << "imgExt=\"jpg\"" << os.widen('\n')
      << "ncams=" << map_viewIdToContiguous.size() << os.widen('\n')
      << "width=" << w_h_image_size.first << os.widen('\n')
      << "height=" << w_h_image_size.second << os.widen('\n')
      << "scale=2" << os.widen('\n')
      << "workDirName=\"_tmp_fast\"" << os.widen('\n')
      << "doPrepareData=FALSE" << os.widen('\n')
      << "doPrematchSifts=FALSE" << os.widen('\n')
      << "doPlaneSweepingSGM=FALSE"  << os.widen('\n')
      << "doFuse=FALSE" << os.widen('\n')
      << os.widen('\n')
      << "[uvatlas]" << os.widen('\n')
      << "texSide=1024" << os.widen('\n')
      << "scale=1" << os.widen('\n')
      << os.widen('\n')
      << "[delanuaycut]" << os.widen('\n')
      << "saveMeshTextured=FALSE" << os.widen('\n')
      << os.widen('\n')
      << "[hallucinationsFiltering]" << os.widen('\n')
      << "useSkyPrior=FALSE" << os.widen('\n')
      << "doLeaveLargestFullSegmentOnly=FALSE" << os.widen('\n')
      << "doRemoveHugeTriangles=TRUE" << os.widen('\n')
      << os.widen('\n')
      << "[largeScale]" << os.widen('\n')
      << "doGenerateAndReconstructSpaceMaxPts=TRUE" << os.widen('\n')
      << "doGenerateSpace=TRUE" << os.widen('\n')
      << "planMaxPts=3000000" << os.widen('\n')
      << "doComputeDEMandOrtoPhoto=FALSE" << os.widen('\n')
      << "doGenerateVideoFrames=FALSE" << os.widen('\n')
      << os.widen('\n')
      << "[meshEnergyOpt]" << os.widen('\n')
      << "doOptimizeOrSmoothMesh=FALSE" << os.widen('\n')
      << os.widen('\n')
      << os.widen('\n')
      << "#EOF" << os.widen('\n')
      << os.widen('\n')
      << os.widen('\n');

    std::ofstream file2(
      stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory),
      "02_mvs_limitedScale" ,"ini").c_str());
    file2 << os.str();
    file2.close();
  }
  return bOk;
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfM_Data_Filename;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch (const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
      << "[-o|--outdir] path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (exportToCMPMVSFormat(sfm_data, stlplus::folder_append_separator(sOutDir) + "CMPMVS"))
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
