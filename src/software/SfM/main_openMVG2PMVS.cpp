
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/SfMViewer/document.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iterator>
#include <iomanip>

bool exportToPMVSFormat(
  const Document & doc,
  const std::string & sOutDirectory,  //Output PMVS files directory
  const std::string & sImagePath  // The images path
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
    // Export data :
    //Camera

    size_t count = 0;
    for (std::map<size_t, PinholeCamera>::const_iterator iter = doc._map_camera.begin();
      iter != doc._map_camera.end(); ++iter, ++count)
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
      iter != doc._map_camera.end();  ++iter, ++count)
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
    os << "level 1" << "\n"
     << "csize 2" << "\n"
     << "threshold 0.7" << "\n"
     << "wsize 7" << "\n"
     << "minImageNum 3" << "\n"
     << "CPU 8" << "\n"
     << "setEdge 0" << "\n"
     << "useBound 0" << "\n"
     << "useVisData 0" << "\n"
     << "sequence -1" << "\n"
     << "timages -1 0 " << doc._map_camera.size() << "\n"
     << "oimages 0" << "\n"; // ?

    std::ofstream file(stlplus::create_filespec(sOutDirectory, "pmvs_options", "txt").c_str());
    file << os.str();
    file.close();
  }
  return bOk;
}

int main(int argc, char *argv[]) {

  Document m_doc;

  // OpenProject
  std::cout << "Type the path to the project :" <<std::endl;
  std::string spath;
  std::cin >> spath;

  std::cout << "\n Open the directory : \n" << spath << std::endl;

  if (m_doc.load(spath))
  {
    exportToPMVSFormat(m_doc,
      stlplus::folder_append_separator(spath) + "pmvs",
      stlplus::folder_append_separator(spath) + "images");

    return( EXIT_SUCCESS );
  }

  // Exit program
  return( EXIT_FAILURE );
}
