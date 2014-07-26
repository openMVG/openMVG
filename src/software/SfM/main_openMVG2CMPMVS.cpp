
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/SfMViewer/document.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>
#include <iomanip>

bool exportToCMPMVSFormat(
  const Document & doc,
  const std::string & sOutDirectory,  //Output CMPMVS files directory
  const std::string & sImagePath  // The images path
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
    // Export data :
    //Camera

    size_t count = 1;
    for (std::map<size_t, PinholeCamera>::const_iterator iter = doc._map_camera.begin();
      iter != doc._map_camera.end(); ++iter, ++count)
    {
      const Mat34 & PMat = iter->second._P;
      std::ostringstream os;
      os << std::setw(5) << std::setfill('0') << count << "_P";
      std::ofstream file(
        stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory),
        os.str() ,"txt").c_str());
      file << "CONTOUR\n"
        << PMat.row(0) <<"\n"<< PMat.row(1) <<"\n"<< PMat.row(2) << std::endl;
      file.close();
    }

    // Image
    count = 1;
	  int w,h; // Image size (suppose they are all the same)
    Image<RGBColor> image;
    for (std::map<size_t, PinholeCamera>::const_iterator iter = doc._map_camera.begin();
      iter != doc._map_camera.end();  ++iter, ++count)
    {
      size_t imageIndex = iter->first;
      const std::string & sImageName = doc._vec_imageNames[imageIndex];
      std::ostringstream os;
      os << std::setw(5) << std::setfill('0') << count;
      ReadImage( stlplus::create_filespec( sImagePath, sImageName).c_str(), &image );
      w = image.Width();
      h = image.Height();
      std::string sCompleteImageName = stlplus::create_filespec(
        stlplus::folder_append_separator(sOutDirectory), os.str(),"jpg");
      WriteImage( sCompleteImageName.c_str(), image);
    }

    // Write the mvs_firstRun script
    std::ostringstream os;
    os << "[global]" << std::endl
    << "dirName=\"" << stlplus::folder_append_separator(sOutDirectory) <<"\"" << std::endl
    << "prefix=\"\"" << std::endl
    << "imgExt=\"jpg\"" << std::endl
    << "ncams=" << doc._map_camera.size() << std::endl
    << "width=" << w << std::endl
    << "height=" << h << std::endl
    << "scale=2" << std::endl
    << "workDirName=\"_tmp_fast\"" << std::endl
    << "doPrepareData=TRUE" << std::endl
    << "doPrematchSifts=TRUE" << std::endl
    << "doPlaneSweepingSGM=TRUE"  << std::endl
    << "doFuse=TRUE" << std::endl
    << "nTimesSimplify=10" << std::endl
    << std::endl
    << "[prematching]" << std::endl
    << "minAngle=3.0" << std::endl
    << std::endl
    << "[grow]" << std::endl
    << "minNumOfConsistentCams=6" << std::endl
    << std::endl
    << "[filter]" << std::endl
    << "minNumOfConsistentCams=2" << std::endl
    << std::endl
    << "#do not erase empy lines after this comment otherwise it will crash ... bug" << std::endl
    << std::endl
    << std::endl;

    std::ofstream file(
	    stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory),
	    "01_mvs_firstRun" ,"ini").c_str());
    file << os.str();
    file.close();

    // limitedScale
    os.str("");
    os << "[global]" << std::endl
    << "dirName=\"" << stlplus::folder_append_separator(sOutDirectory) <<"\"" << std::endl
    << "prefix=\"\"" << std::endl
    << "imgExt=\"jpg\"" << std::endl
    << "ncams=" << doc._map_camera.size() << std::endl
    << "width=" << w << std::endl
    << "height=" << h << std::endl
    << "scale=2" << std::endl
    << "workDirName=\"_tmp_fast\"" << std::endl
    << "doPrepareData=FALSE" << std::endl
    << "doPrematchSifts=FALSE" << std::endl
    << "doPlaneSweepingSGM=FALSE"  << std::endl
    << "doFuse=FALSE" << std::endl
    << std::endl
    << "[uvatlas]" << std::endl
    << "texSide=1024" << std::endl
    << "scale=1" << std::endl
    << std::endl
    << "[delanuaycut]" << std::endl
    << "saveMeshTextured=FALSE" << std::endl
    << std::endl
    << "[hallucinationsFiltering]" << std::endl
    << "useSkyPrior=FALSE" << std::endl
    << "doLeaveLargestFullSegmentOnly=FALSE" << std::endl
    << "doRemoveHugeTriangles=TRUE" << std::endl
    << std::endl
    << "[largeScale]" << std::endl
    << "doGenerateAndReconstructSpaceMaxPts=TRUE" << std::endl
    << "doGenerateSpace=TRUE" << std::endl
    << "planMaxPts=3000000" << std::endl
    << "doComputeDEMandOrtoPhoto=FALSE" << std::endl
    << "doGenerateVideoFrames=FALSE" << std::endl
    << std::endl
    << "[meshEnergyOpt]" << std::endl
    << "doOptimizeOrSmoothMesh=FALSE" << std::endl
    << std::endl
    << std::endl
    << "#EOF" << std::endl
    << std::endl
    << std::endl;

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
  std::string sSfMDir;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfMDir, "sfmdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdir path, the SfM_output path]\n"
      << "[-o|--outdir path]\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );
  
  Document m_doc;
  if (m_doc.load(sSfMDir))
  {
    exportToCMPMVSFormat(m_doc,
      stlplus::folder_append_separator(sOutDir) + "CMPMVS",
      stlplus::folder_append_separator(sSfMDir) + "images");

    return( EXIT_SUCCESS );
  }

  // Exit program
  return( EXIT_FAILURE );
}
