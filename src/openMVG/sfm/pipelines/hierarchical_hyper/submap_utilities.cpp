// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_io_cereal.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_utilities.hpp"

#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>
#include <cereal/archives/json.hpp>

#include <iostream>
#include <fstream>

namespace openMVG{
namespace sfm{

bool SaveHsfmSubmap(
    const HsfmSubmap & submap,
    const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  try
  {
    cereal::JSONOutputArchive archive(stream);
    archive(cereal::make_nvp("submap", submap));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }

  stream.close();
  return true;
}

bool LoadHsfmSubmap(
    HsfmSubmap & submap,
    const std::string & filename)
{
  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  try
  {
    cereal::JSONInputArchive archive(stream);
    archive(cereal::make_nvp("submap", submap));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }

  stream.close();
  return true;
}

bool SaveSubmaps(
    const HsfmSubmaps & submaps,
    const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  try
  {
    cereal::JSONOutputArchive archive(stream);
    archive(cereal::make_nvp("submaps", submaps));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }

  stream.close();
  return true;
}


bool LoadSubmaps(HsfmSubmaps &submaps,
    const std::string & filename)
{
  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  try
  {
    cereal::JSONInputArchive archive(stream);
    archive(cereal::make_nvp("submaps", submaps));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }

  stream.close();
  return true;
}

/**
 * @brief Finds the sibling submap id
 * @param submaps
 * @param submap id
 * @return sibling id
 */
IndexT getSiblingSubmapId(const HsfmSubmaps & submaps, const IndexT submap_id)
{
  const IndexT parent_id = submaps.at(submap_id).parent_id;
  std::pair<IndexT,IndexT> siblings_submaps_ids = submaps.at(parent_id).children_submaps;
  if (siblings_submaps_ids.first == submap_id)
    return siblings_submaps_ids.second;
  return siblings_submaps_ids.first;
}

/**
 * @brief checks a submap quality
 * @param submaps A maps of all submaps
 * @param submap id of the current submap
 * @return a diagnostic of the submap
 * @note checks number of poses/number of views and how many separator tracks are reconstructed,
 * as well as how many of these tracks have a common reconstruction in the sibling submap
 */
eSubmapDiagnostic CheckSingleSubmap(const HsfmSubmaps & submaps, const IndexT & submap_id)
{
  // diagnostic
  eSubmapDiagnostic Diagnostic = OK;

  // *******************************************
  // Start by checking reconstruction, i.e.
  // if there are enough reconstructed poses.
  // *******************************************
  
  const HsfmSubmap & submap = submaps.at(submap_id);
  const SfM_Data & sfm_data = submap.sfm_data;
  const int n_poses = sfm_data.GetPoses().size();
  const int n_views = sfm_data.GetViews().size();

  std::cout << "\nSUBMAP # " << submap_id << " diagnostic : " << std::endl;
  std::cout << "number of poses : " << n_poses << std::endl;
  std::cout << "number of views : " << n_views << std::endl;
  std::cout << "ratio : " << 100*(n_poses/(double)n_views) << " %" << std::endl;
  std::cout << " === \n";

  if (n_poses <= 2)
  {
    Diagnostic = NO_RECONSTRUCTION;
    return Diagnostic;
  }
  else if (n_poses/(double)n_views < 0.3) // TODO : change this value to a parameter ! 
  {
    Diagnostic = BAD_RECONSTRUCTION;
    return Diagnostic;
  }

  // *******************************************
  // if reconstruction is ok, check that we 
  // have enough reconstructed separators 
  // *******************************************
  
  const IndexT parent_id = submap.parent_id;
  const std::set<IndexT> & separator_tracks = submaps.at(parent_id).separator;
  const int n_total_separators = separator_tracks.size();
  std::set<IndexT> reconstructed_separators;

  // count reconstructed separators
  for (const auto & landmark : sfm_data.GetLandmarks())
  {
    const IndexT & landmark_id = landmark.first;
    if (separator_tracks.count(landmark_id) != 0)
    {
      reconstructed_separators.insert(landmark_id);
    }
  }
  const int n_reconstructed_separators = reconstructed_separators.size();
  std::cout << "total number of separators : " << n_total_separators << std::endl;
  std::cout << "number of reconstructed separators : " << n_reconstructed_separators << std::endl;
  std::cout << "ratio : " << 100*(n_reconstructed_separators/(double)n_total_separators) << " %" << std::endl;
  std::cout << std::endl;

  if (n_reconstructed_separators == 0)
  {
    Diagnostic = NO_RECONSTRUCTED_SEPARATORS;
    return Diagnostic;
  }
  // TODO : add a NOT_ENOUGH_RECONSTRUCTED_SEPARATORS thing ? whats the threshold there ? :/

  // *******************************************
  // now check that siblings submaps have 
  // enough common reconstructed separators for
  // a clean merging.
  // *******************************************

  const IndexT sister_submap_id = getSiblingSubmapId(submaps, submap_id);

  const SfM_Data & sister_sfm_data = submaps.at(sister_submap_id).sfm_data;
  std::set<size_t> common_reconstructed_separators;
  for (const auto & landmark : sister_sfm_data.GetLandmarks())
  {
    const IndexT & landmark_id = landmark.first;
    if (reconstructed_separators.count(landmark_id) != 0)
    {
      common_reconstructed_separators.insert(landmark_id);
    }
  }
  
  const int n_common_reconstructed_separators = common_reconstructed_separators.size();
  std::cout << "number of common reconstructed separators with sister submap : " 
    << n_common_reconstructed_separators << std::endl;
  std::cout << "total number of common separators with sister submap :         " 
    << n_total_separators << std::endl;
  std::cout << "part of total n of separators :                                " 
    << 100*(n_common_reconstructed_separators/(double)n_total_separators) << "%" << std::endl;
  std::cout << "part of reconstructed n of separators :                        " 
    << 100*(n_common_reconstructed_separators/(double)n_reconstructed_separators) << "%" << std::endl;
  std::cout << std::endl;

  if (n_common_reconstructed_separators == 0)
  {
    Diagnostic = NO_COMMON_RECONSTRUCTED_SEPARATORS;
    return Diagnostic;
  }

  return Diagnostic;
}

/**
 * @brief check all submaps
 * @param submaps
 * @return A map of diagnostics for each of the submaps inputed
 */
std::map<IndexT, eSubmapDiagnostic> CheckSubmaps(const HsfmSubmaps &submaps)
{
  std::map<IndexT, eSubmapDiagnostic> diagnostics;
  for ( const auto & smap : submaps)
  {
    const IndexT & smap_id = smap.first;
    const HsfmSubmap & submap = smap.second;

    eSubmapDiagnostic diagnostic = CheckSingleSubmap(submaps, smap_id);
    diagnostics[smap_id] = diagnostic;
  }

  return diagnostics;
}

}//namespace sfm
}//namespace openMVG
