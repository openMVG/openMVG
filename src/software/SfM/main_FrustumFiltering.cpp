// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/frustum.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_filters_frustum.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

/// Build a list of pair that share visibility content from the SfM_Data structure
Pair_Set BuildPairsFromStructureObservations(const SfM_Data & sfm_data)
{
  Pair_Set pairs;

  for (Landmarks::const_iterator itL = sfm_data.GetLandmarks().begin();
    itL != sfm_data.GetLandmarks().end(); ++itL)
  {
    const Landmark & landmark = itL->second;
    for (Observations::const_iterator iterI = landmark.obs.begin();
      iterI != landmark.obs.end(); ++iterI)
    {
      const IndexT id_viewI = iterI->first;
      Observations::const_iterator iterJ = landmark.obs.begin();
      std::advance(iterJ, 1);
      for (; iterJ != landmark.obs.end(); ++iterJ)
      {
        const IndexT id_viewJ = iterJ->first;
        pairs.insert( std::make_pair(id_viewI,id_viewJ));
      }
    }
  }
  return pairs;
}

/// Build a list of pair from the camera frusta intersections
Pair_Set BuildPairsFromFrustumsIntersections(
  const SfM_Data & sfm_data,
  const double z_near = -1., // default near plane
  const double z_far = -1.,  // default far plane
  const std::string & sOutDirectory = "") // output directory to save frustums as PLY
{
  const Frustum_Filter frustum_filter(sfm_data, z_near, z_far);
  if (!sOutDirectory.empty())
    frustum_filter.export_Ply(stlplus::create_filespec(sOutDirectory, "frustums.ply"));
  return frustum_filter.getFrustumIntersectionPairs();
}

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Frustum filtering:\n"
    << "-----------------------------------------------------------\n"
    << "Compute camera cones that share some putative visual content.\n"
    << "------------------------------------------------------------"
    << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutFile;
  double z_near = -1.;
  double z_far = -1.;

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('o', sOutFile, "output_file") );
  cmd.add( make_option('n', z_near, "z_near") );
  cmd.add( make_option('f', z_far, "z_far") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-o|--output_file] filename of the output pair file\n"
    << "[-n|--z_near] 'optional' distance of the near camera plane\n"
    << "[-f|--z_far] 'optional' distance of the far camera plane\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Assert that we can create the output directory
  if (!stlplus::folder_exists( stlplus::folder_part(sOutFile)) &&
      !stlplus::folder_create( stlplus::folder_part(sOutFile)))
      return EXIT_FAILURE;

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  openMVG::system::Timer timer;

  const Pair_Set pairs = BuildPairsFromFrustumsIntersections(sfm_data, z_near, z_far, stlplus::folder_part(sOutFile));
  /*const Pair_Set pairs = BuildPairsFromStructureObservations(sfm_data); */

  std::cout << "#pairs: " << pairs.size() << std::endl;
  std::cout << std::endl << " Pair filtering took (s): " << timer.elapsed() << std::endl;

  // export pairs on disk
  if (savePairs(sOutFile, pairs))
  {
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
