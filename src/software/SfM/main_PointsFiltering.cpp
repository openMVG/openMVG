// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <memory>
#include <string>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

/**
* @brief Implement a statistical Structure filter that remove 3D points that have:
* - a depth that is too large (threshold computed as factor * median ~= X84)
* @param sfm_data The sfm scene to filter (inplace filtering)
* @param k_factor The factor applied to the median depth per view
* @param k_min_point_per_pose Keep only poses that have at least this amount of points
* @param k_min_track_length Keep only tracks that have at least this length
* @return The min_median_value observed for all the view
*/
double DepthCleaning
(
  SfM_Data & sfm_data,
  const double k_factor =  5.2,    //5.2 * median ~= X84,
  const IndexT k_min_point_per_pose = 12,  // 6 min
  const IndexT k_min_track_length = 2      // 2 min
)
{
  using DepthAccumulatorT = std::vector<double>;
  std::map<IndexT, DepthAccumulatorT > map_depth_accumulator;

  // For each landmark accumulate the camera/point depth info for each view
  for (const auto & landmark_it : sfm_data.structure)
  {
    const Observations & obs = landmark_it.second.obs;
    for (const auto & obs_it : obs)
    {
      const View * view = sfm_data.views.at(obs_it.first).get();
      if (sfm_data.IsPoseAndIntrinsicDefined(view))
      {
        const Pose3 pose = sfm_data.GetPoseOrDie(view);
        const double depth = pose.depth(landmark_it.second.X);
        if (depth > 0)
        {
          map_depth_accumulator[view->id_view].push_back(depth);
        }
      }
    }
  }

  double min_median_value = std::numeric_limits<double>::max();
  std::map<IndexT, double > map_median_depth;
  for (const auto & iter : sfm_data.GetViews())
  {
    const View * v = iter.second.get();
    const IndexT view_id = v->id_view;
    if (map_depth_accumulator.count(view_id) == 0)
      continue;
    // Compute median from the depth distribution
    const auto & acc = map_depth_accumulator.at(view_id);
    double min, max, mean, median;
    if (minMaxMeanMedian(acc.begin(), acc.end(), min, max, mean, median))
    {

      min_median_value = std::min(min_median_value, median);
      // Compute depth threshold for each view: factor * medianDepth
      map_median_depth[view_id] = k_factor * median;
    }
  }
  map_depth_accumulator.clear();

  // Delete invalid observations
  size_t cpt = 0;
  for (auto & landmark_it : sfm_data.structure)
  {
    Observations obs;
    for (auto & obs_it : landmark_it.second.obs)
    {
      const View * view = sfm_data.views.at(obs_it.first).get();
      if (sfm_data.IsPoseAndIntrinsicDefined(view))
      {
        const Pose3 pose = sfm_data.GetPoseOrDie(view);
        const double depth = pose.depth(landmark_it.second.X);
        if ( depth > 0
            && map_median_depth.count(view->id_view)
            && depth < map_median_depth[view->id_view])
          obs.insert(obs_it);
        else
          ++cpt;
      }
    }
    landmark_it.second.obs.swap(obs);
  }
  std::cout << "#point depth filter: " << cpt << " measurements removed" <<std::endl;

  // Remove orphans
  eraseUnstablePosesAndObservations(sfm_data, k_min_point_per_pose, k_min_track_length);

  return min_median_value;
}

} // namespace sfm
} // namespace openMVG

using namespace openMVG;
using namespace openMVG::sfm;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Filter Structure based on statistics computed per view    :\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sfm_data_filename;
  std::string sfm_data_filename_out;
  double factor = 5.2; // 5.2 * median ~= X84

  cmd.add( make_option('i', sfm_data_filename, "input_file") );
  cmd.add( make_option('o', sfm_data_filename_out, "output_file") );
  cmd.add( make_option('f', factor, "factor") );

  try
  {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string& s)
  {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--input_file] path to a SfM_Data scene\n"
      << "[-o|--output_file] path where the filtered SfM_data scene will be saved\n"
      << "\n[Optional]\n"
      << "[-f|--factor] factor apply on the median depth per view for thresholding\n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sfm_data_filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sfm_data_filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Structure filtering
  //---------------------------------------

  const IndexT k_min_point_per_pose = 12;
  const IndexT k_min_track_length = 2;
  const double k_min_median_depth = DepthCleaning(
    sfm_data,
    factor,
    k_min_point_per_pose,
    k_min_track_length
  );
  std::cout << "MIN MEDIAN DEPTH VALUE = " << k_min_median_depth << std::endl;

  if (!Save(sfm_data, sfm_data_filename_out, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The output SfM_Data file \""<< sfm_data_filename_out << "\" cannot be saved." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
