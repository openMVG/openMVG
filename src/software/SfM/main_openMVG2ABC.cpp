// Copyright (c) 2015 cpichard.
//
// // This Source Code Form is subject to the terms of the Mozilla Public
// // License, v. 2.0. If a copy of the MPL was not distributed with this
// // file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <cstdlib>

// Command line
#include "third_party/cmdLine/cmdLine.h"

// OpenMVG
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/geometry/pose3.hpp"
#include <openMVG/dataio/AlembicExporter.hpp>
using namespace openMVG;

// Alembic
//#include <Alembic/AbcGeom/All.h>
//#include <Alembic/AbcCoreHDF5/All.h>
using namespace Alembic::Abc;
namespace AbcG = Alembic::AbcGeom;
using namespace AbcG;

int main(int argc, char **argv)
{

  // Get arguments
  CmdLine cmdLine;
  std::string sfmDataFilename;
  std::string sOutAlembic = "";

  cmdLine.add(make_option('i', sfmDataFilename, "sfmdata"));
  cmdLine.add(make_option('o', sOutAlembic, "outfile"));

  try
  {
    if(argc < 4) throw std::string("Invalid number of parameters in the command line.");
    cmdLine.process(argc, argv);
  }
  catch(const std::string &s)
  {
    std::cout << "openMVG to alembic\n";
    std::cout << "Usage: " << argv[0] << '\n'
            << "[-i|--sfmdata filename, the SfM_Data file to convert]\n"
            << "[-o|--outfile path]\n"
            << std::endl;
    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  //
  sfm::SfM_Data sfm_data;
  if(!Load(sfm_data, sfmDataFilename, sfm::ESfM_Data(sfm::ESfM_Data::ALL)))
  {
    std::cout << "Error: The input project file \""
            << sfmDataFilename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Open Alembic archive
  // TODO: decide wether we want to keep HDF5 or switch to Ogawa 
  dataio::AlembicExporter exporter(sOutAlembic);

  // Export points
  exporter.addPoints(sfm_data.GetLandmarks());
 
  // Export cameras
  exporter.addCameras(sfm_data);

  return EXIT_SUCCESS;
}
