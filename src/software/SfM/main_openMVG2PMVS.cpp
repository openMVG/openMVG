
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/SfMViewer/document.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>
#include <iomanip>

bool exportToPMVSFormat(
  const Document & doc,
  const std::string & sOutDirectory,  //Output PMVS files directory
  const std::string & sImagePath,  // The images path
  const int resolution,
  const int CPU
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
    C_Progress_display my_progress_bar( doc._map_camera.size()*2 );
    // Export data :
    //Camera

    size_t count = 0;
    for (std::map<size_t, PinholeCamera>::const_iterator iter = doc._map_camera.begin();
      iter != doc._map_camera.end(); ++iter, ++count, ++my_progress_bar)
    {
      const Mat34 & PMat = iter->second._P;
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << count;
      std::ofstream file(
        stlplus::create_filespec(stlplus::folder_append_separator(sOutDirectory) + "txt",
        os.str() ,"txt").c_str());
      file << "CONTOUR\n"
        << PMat.row(0) <<"\n"<< PMat.row(1) <<"\n"<< PMat.row(2) << std::endl;
      file.close();
    }

    // Image
    count = 0;
    Image<RGBColor> image;
    for (std::map<size_t, PinholeCamera>::const_iterator iter = doc._map_camera.begin();
      iter != doc._map_camera.end();  ++iter, ++count, ++my_progress_bar)
    {
      size_t imageIndex = iter->first;
      const std::string & sImageName = doc._vec_imageNames[imageIndex];
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << count;
      ReadImage( stlplus::create_filespec( sImagePath, sImageName).c_str(), &image );
      std::string sCompleteImageName = stlplus::create_filespec(
        stlplus::folder_append_separator(sOutDirectory) + "visualize", os.str(),"jpg");
      WriteImage( sCompleteImageName.c_str(), image);
    }

    //pmvs_options.txt
    std::ostringstream os;
    os << "level " << resolution << "\n"
     << "csize 2" << "\n"
     << "threshold 0.7" << "\n"
     << "wsize 7" << "\n"
     << "minImageNum 3" << "\n"
     << "CPU " << CPU << "\n"
     << "setEdge 0" << "\n"
     << "useBound 0" << "\n"
     << "useVisData 0" << "\n"
     << "sequence -1" << "\n"
     << "maxAngle 10" << "\n"
     << "quad 2.0" << "\n"
     << "timages -1 0 " << doc._map_camera.size() << "\n"
     << "oimages 0" << "\n"; // ?

    std::ofstream file(stlplus::create_filespec(sOutDirectory, "pmvs_options", "txt").c_str());
    file << os.str();
    file.close();
  }
  return bOk;
}

int main(int argc, char *argv[]) {

  CmdLine cmd;
  std::string sSfMDir;
  std::string sOutDir = "";
  int resolution = 1;
  int CPU = 8;

  cmd.add( make_option('i', sSfMDir, "sfmdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', resolution, "resolution") );
  cmd.add( make_option('c', CPU, "CPU") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--sfmdir path, the SfM_output path] "
      << "[-o|--outdir path] "
      << "[-r|--resolution: divide image coefficient] "
      << "[-c|--nb core] "
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
    exportToPMVSFormat(m_doc,
      stlplus::folder_append_separator(sOutDir) + "PMVS",
      stlplus::folder_append_separator(sSfMDir) + "images",
      resolution,
      CPU );

    return( EXIT_SUCCESS );
  }

  // Exit program
  return( EXIT_FAILURE );
}
