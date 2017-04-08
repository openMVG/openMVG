// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre Moulon, Srivathsan Murali, Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/sfm/sfm.hpp"
#include "openMVG/system/timer.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "domSetLibrary/domset.h"
#include "domSetLibrary/types.h"

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

/**
* @brief convert openMVG sfm_data to domset data
* @param sfm_data openMVG dataset
* @param[out] cameras list of intrinsics
* @param[out] view list of views
* @param[out] points list of points
* @param[out] map_view map between views added and the global index
*  We need this map to handle the fact that OpenMVG view Id can be non
*   contiguous:
*   - in case of missing pose or intrinsic for some view
* @retval true if success
* @retval false if failure
*/
bool domsetImporter(
  const SfM_Data &sfm_data,
  std::vector<nomoko::Camera> &cameras,
  std::vector<nomoko::View> &views,
  std::vector<nomoko::Point> &points,
  std::map<openMVG::IndexT, uint32_t> &map_view )
{

  // Convert OpenMVG data to domset library data
  openMVG::system::Timer loadDataTimer;

  // adding views
  for ( const auto &view : sfm_data.GetViews() )
  {
    if ( sfm_data.IsPoseAndIntrinsicDefined( view.second.get() ) )
    {
      map_view[ view.first ] = views.size(); // need to make global

      const openMVG::geometry::Pose3 poseMVG( sfm_data.GetPoseOrDie( view.second.get() ) );
      nomoko::View v;
      v.rot   = poseMVG.rotation().cast<float>();
      v.trans = poseMVG.center().transpose().cast<float>();
      views.push_back( v );
    }
  }

  // adding landmarks
  for ( const auto &it_landmark : sfm_data.GetLandmarks())
  {
    const Landmark &landmark = it_landmark.second;
    const Observations &obs  = landmark.obs;
    // X, color, obsCount
    std::vector<size_t> vIds;
    for ( const auto &it_obs : obs)
    {
      vIds.push_back( map_view[ it_obs.first ] );
    }

    nomoko::Point p;
    p.pos      = landmark.X.transpose().cast<float>();
    p.viewList = vIds;
    points.push_back( p );
  }

  std::cout << std::endl
            << "Number of views  = " << views.size() << std::endl
            << "Number of points = " << points.size() << std::endl
            << "Loading data took (s): "
            << loadDataTimer.elapsed() << std::endl;
  return true;
}

/**
* @brief Export a sfm_data file using a subset of the view of a given sfm_data
* @param sfm_data The whole data set
* @param outFilename Output file name
* @param cluster List of view to consider
* @retval true if success
* @retval false if failure
*/
bool exportData( const SfM_Data &sfm_data,
                 const std::string &outFilename,
                 const std::set<size_t> &cluster )
{
  SfM_Data cl_sfm_data;
  cl_sfm_data.s_root_path = sfm_data.s_root_path;

  // Copy the view (only the requested ones)
  for ( const auto view : sfm_data.GetViews() )
  {
    const bool inCluster = cluster.find(view.first) != cluster.end();
    if ( inCluster && sfm_data.IsPoseAndIntrinsicDefined( view.second.get() ) )
    {
      cl_sfm_data.poses[ view.first ] = sfm_data.GetPoseOrDie( view.second.get() );
      cl_sfm_data.views[ view.first ] = view.second;

      const auto intrinsic  = sfm_data.GetIntrinsics().at( view.second.get()->id_intrinsic );
      if (cl_sfm_data.intrinsics.count(view.second.get()->id_intrinsic) == 0)
        cl_sfm_data.intrinsics[ view.second.get()->id_intrinsic ] = intrinsic;
    }
  }

  // Copy observations that have relation with the considered view
  for ( const auto &it_landmark : sfm_data.GetLandmarks() )
  {
    const Landmark &landmark = it_landmark.second;
    Observations obs;
    for ( const auto &observation : landmark.obs )
    {
      const auto &it = cl_sfm_data.views.find( observation.first );
      if ( it != cl_sfm_data.views.end() )
      {
        obs[ observation.first ] = observation.second;
      }
    }
    // Landmark observed in less than 3 view are ignored
    if ( obs.size() < 2 )
      continue;
    cl_sfm_data.structure[ it_landmark.first ].X   = landmark.X;
    cl_sfm_data.structure[ it_landmark.first ].obs = obs;
  }

  return Save( cl_sfm_data, outFilename, ESfM_Data( ALL ) );
}

int main( int argc, char **argv )
{
  using namespace std;
  std::cout << "Dominant Set Clustering" << std::endl
            << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename     = "";
  std::string sOutDir                = "";
  unsigned int clusterSizeLowerBound = 20;
  unsigned int clusterSizeUpperBound = 30;
  float voxelGridSize                = 10.0f;

  cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  cmd.add( make_option( 'o', sOutDir, "outdir" ) );
  cmd.add( make_option( 'l', clusterSizeLowerBound, "cluster_size_lower_bound" ) );
  cmd.add( make_option( 'u', clusterSizeUpperBound, "cluster_size_upper_bound" ) );
  cmd.add( make_option( 'v', voxelGridSize, "voxel_grid_size" ) );

  try
  {
    if ( argc == 1 )
      throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch ( const std::string &s )
  {
    std::cerr << "Usage: " << argv[ 0 ] << "\n"
              << "[-i|--input_file] path to a SfM_Data scene\n"
              << "[-o|--outdir path] path to output directory\n"
              << "[-l|--cluster_size_lower_bound] lower bound to cluster size\n"
              << "[-u|--cluster_size_upper_bound] upper bound to cluster size\n"
              << "[-v|--voxel_grid_size] voxel grid size\n"
              << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Params: " << argv[ 0 ]  << std::endl
            << "[Input file]      = "   << sSfM_Data_Filename << std::endl
            << "[Outdir path]     = "   << sOutDir << std::endl
            << "[Cluster size:"         << std::endl
            << "    Lower bound   = "   << clusterSizeLowerBound << std::endl
            << "    Upper bound]   = "  << clusterSizeUpperBound << std::endl
            << "[Voxel grid size]  = "  << voxelGridSize << std::endl;

  if ( sSfM_Data_Filename.empty() )
  {
    std::cerr << "\nIt is an invalid file input" << std::endl;
    return EXIT_FAILURE;
  }

  // Prepare output folder
  if ( !stlplus::folder_exists( sOutDir ) )
    if ( !stlplus::folder_create( sOutDir ))
    {
      std::cerr << "Cannot create: " << sOutDir << std::endl;
      return EXIT_FAILURE;
    }

  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfM_Data_Filename, ESfM_Data( ALL ) ) )
  {
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfM_Data_Filename << "\" can't be read."
              << std::endl;
    return EXIT_FAILURE;
  }

  // loading data
  std::vector<nomoko::Camera> cameras; // stores the various camera intrinsic parameters
  std::vector<nomoko::View> views;     // stores the poses for each view
  std::vector<nomoko::Point> points;   // 3d point positions

  std::map<openMVG::IndexT, uint32_t> origViewMap; // need to keep track of original views ids
  if ( !domsetImporter( sfm_data, cameras, views, points, origViewMap ) )
  {
    std::cerr << "Error: can't import data" << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // clustering views process
  //---------------------------------------
  openMVG::system::Timer clusteringTimer;

  nomoko::Domset domset( points, views, cameras, voxelGridSize );
  domset.clusterViews( clusterSizeLowerBound, clusterSizeUpperBound );

  std::cout << "Clustering view took (s): "
            << clusteringTimer.elapsed() << std::endl;

  // export to ply to visualize
  const std::string viewOut = sOutDir + "/views.ply";
  domset.exportToPLY( viewOut );

  // Retrieve the cluster and export them
  std::vector<std::set<size_t>> finalClusters;
  {
    const std::vector<std::vector<size_t>> clusters = domset.getClusters();

    // Remap the camera index from contiguous to original view Id
    std::map<openMVG::IndexT, uint32_t> origViewMap_reverse;
    for (const auto& it : origViewMap) {
      origViewMap_reverse.insert(std::make_pair(it.second, it.first));
    }
    // For every cluster, remap the view Id
    for ( const auto &cl : clusters )
    {
      std::set<size_t> newCl;
      for ( const auto vId : cl )
        newCl.insert( origViewMap_reverse[ vId ] );
      finalClusters.emplace_back(newCl);
    }
  }

  const size_t numClusters = finalClusters.size();
  std::cout << "Number of clusters = " << numClusters << std::endl;

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for
#endif
  for ( int i = 0; i < numClusters; ++i )
  {
    std::stringstream filename;
    filename << sOutDir << "/sfm_data";
    filename << std::setw( 4 ) << std::setfill( '0' ) << i;
    filename << ".bin";

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
    {
      std::stringstream ss;
      ss << "Writing cluster to " << filename.str() << std::endl;
      std::cout << ss.str();
    }

    if ( !exportData( sfm_data, filename.str(), finalClusters[ i ] ) )
    {
      std::stringstream str;
      str << "Could not write cluster : " << filename.str() << std::endl;
      std::cerr << str.str();
    }
  }
  return EXIT_SUCCESS;
}
