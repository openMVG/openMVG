
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_IO_H
#define OPENMVG_SFM_IO_H

#include <vector>
#include <string>
#include <fstream>
#include <iterator>

namespace openMVG{
namespace SfMIO{

// Load an image file list
// One basename per line.
bool loadImageList(std::vector<std::string> & vec_fileNames,
  std::string sFileName = "./lists.txt",
  bool bVerbose = true)
{
  std::ifstream in(sFileName.c_str());
  if(!in.is_open())  {
    std::cerr << std::endl
      << "Impossible to read the specified file." << std::endl;
  }
  std::string sValue;
  while(in>>sValue)
    vec_fileNames.push_back(sValue);

  // DEBUG INFO
  if (bVerbose) {
    std::cout << "\nIMAGE :" << std::endl;
    std::copy(vec_fileNames.begin(), vec_fileNames.end(),
      std::ostream_iterator<std::string>(std::cout, "\n"));
  }
  in.close();
  return !(vec_fileNames.empty());
}

bool loadIntrinsic(const std::string & fileName, Mat3 & K)
{
  // Load the K matrix
  std::ifstream in;
  in.open( fileName.c_str(), std::ifstream::in);
  if(in.is_open())  {
    for (int j=0; j < 3; ++j)
      for (int i=0; i < 3; ++i)
        in >> K(j,i);
  }
  else  {
    std::cerr << std::endl
      << "Invalid input K.txt file" << std::endl;
    return false;
  }
  bool bOk = !in.bad();
  in.close();
  return bOk;
}


} // namespace SfMIO
} // namespace openMVG

#endif // OPENMVG_SFM_INCREMENTAL_ENGINE_H

