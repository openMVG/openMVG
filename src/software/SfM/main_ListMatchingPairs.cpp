// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress_display.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <cstdlib>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;

enum ePairMode
{
  PAIR_MODE_EXHAUSTIVE = 0,
  PAIR_MODE_CONTIGUOUS = 1,
  PAIR_MODE_NEIGHBORHOOD = 2
};

/// Export an adjacency matrix as a SVG file
void AdjacencyMatrixToSVG
(
  const size_t NbImages,
  const Pair_Set & corresponding_indexes,
  const std::string & sOutName
)
{
  using namespace svg;
  if (!corresponding_indexes.empty())
  {
    const float scaleFactor = 5.0f;
    svgDrawer svgStream((NbImages+3)*5, (NbImages+3)*5);
    // List possible pairs
    for (size_t I = 0; I < NbImages; ++I)
    {
      for (size_t J = 0; J < NbImages; ++J)
      {
        // If the pair have matches display a blue boxes at I,J position.
        const auto iterSearch = corresponding_indexes.find(std::make_pair(I,J));
        if (iterSearch != corresponding_indexes.end())
        {
          svgStream.drawSquare(J*scaleFactor, I*scaleFactor, scaleFactor/2.0f,
          svgStyle().fill("blue").noStroke());
        }
      }
    }
    // Display axes with 0 -> NbImages annotation : _|
    std::ostringstream osNbImages;
    osNbImages << NbImages;
    svgStream.drawText((NbImages+1)*scaleFactor, scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine((NbImages+1)*scaleFactor, 2*scaleFactor,
      (NbImages+1)*scaleFactor, (NbImages)*scaleFactor - 2*scaleFactor,
      svgStyle().stroke("black", 1.0));

    svgStream.drawText(scaleFactor, (NbImages+1)*scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages)*scaleFactor - scaleFactor,
      (NbImages+1)*scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine(2*scaleFactor, (NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - 2*scaleFactor, (NbImages+1)*scaleFactor,
      svgStyle().stroke("black", 1.0));

    std::ofstream svgFileStream(sOutName.c_str());
    svgFileStream << svgStream.closeSvgFile().str();
  }
}

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Compute a view pair list file for main_ComputeMatches:\n"
    << " - various pair modes are available to adapt to user dataset\n"
    << "   configuration.\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string s_SfM_Data_filename;
  std::string s_out_file;
  int i_neighbor_count = 5;
  int i_mode(PAIR_MODE_EXHAUSTIVE);

  cmd.add( make_option('i', s_SfM_Data_filename, "input_file") );
  cmd.add( make_option('o', s_out_file, "output_file") );
  cmd.add( make_option('n', i_neighbor_count, "neighbor_count") );
  cmd.add( make_switch('G', "gps_mode"));
  cmd.add( make_switch('V', "video_mode"));
  cmd.add( make_switch('E', "exhaustive_mode"));

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-o|--output_file] the output pairlist file (i.e ./pair_list.txt)\n"
    << "optional:\n"
    << "Matching pair modes [E/V/G]:\n"
    << "\t[-E|--exhaustive_mode] exhaustive mode (default mode)\n"
    << "\t[-V|--video_mode] link views that belongs to contiguous poses ids\n"
    << "\t[-G|--gps_mode] use the pose center priors to link neighbor views\n"
    << "Note: options V & G are linked the following parameter:\n"
    << "\t [-n|--neighbor_count] number of maximum neighbor\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout
    << " You called : " << "\n"
    << argv[0] << "\n"
    << "--input_file " << s_SfM_Data_filename << "\n"
    << "--output_file " << s_out_file << "\n"
    << "Optional parameters:" << "\n"
    << "--exhaustive_mode " << (cmd.used('E') ? "ON" : "OFF") << "\n"
    << "--video_mode " <<  (cmd.used('V') ? "ON" : "OFF") << "\n"
    << "--gps_mode "  << (cmd.used('G') ? "ON" : "OFF") << "\n";
  if (cmd.used('V') || cmd.used('G'))
    std::cout << "--neighbor_count " << i_neighbor_count << std::endl;

  std::cout << std::endl;

  //--
  // Check validity of the input parameters
  //--

  // pair list mode
  if ( int(cmd.used('E')) + int(cmd.used('V')) + int(cmd.used('G')) > 1)
  {
    std::cerr << "You can use only one matching mode." << std::endl;
    return EXIT_FAILURE;
  }
  if (cmd.used('E'))
    i_mode = PAIR_MODE_EXHAUSTIVE;
  else if (cmd.used('V'))
    i_mode = PAIR_MODE_CONTIGUOUS;
  else if (cmd.used('G'))
    i_mode = PAIR_MODE_NEIGHBORHOOD;

  // Input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, s_SfM_Data_filename, ESfM_Data(VIEWS|INTRINSICS)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << s_SfM_Data_filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout
    << "Loaded a sfm_data scene with:\n"
    << " #views: " << sfm_data.GetViews().size() << "\n"
    << std::endl;

  // out file
  if (s_out_file.empty())
  {
    std::cerr << "Invalid output filename." << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(stlplus::folder_part(s_out_file)))
  {
    if (!stlplus::folder_create(stlplus::folder_part(s_out_file)))
    {
      std::cerr << "Cannot create directory for the output file." << std::endl;
      return EXIT_FAILURE;
    }
  }

  //---------------------------------------
  // a. List the view pose as a linear sequence of ids.
  // b. Establish a pose graph according the user chosen mode:
  //    - E => upper diagonal pairs,
  //    - V => list the N closest pose ids,
  //    - G => list the N closest poses XYZ position.
  // c. Convert the pose graph edges to a view graph
  // d. Export the view graph to a file and a SVG adjacency list
  //---------------------------------------
  std::multimap<IndexT, IndexT> pose_id_toViewId;
  std::set<IndexT> set_poses;

  // a. Get nodes of the pose graph as a linear sequence
  for (const auto & viewIter : sfm_data.GetViews())
  {
    const View * v = viewIter.second.get();
    assert (viewIter.first == v->id_view);
    pose_id_toViewId.insert( std::make_pair(v->id_pose, v->id_view) );
    set_poses.insert(v->id_pose);
  }
  const std::vector<IndexT> vec_poses(set_poses.begin(), set_poses.end());


  // b. Create the pose graph pair relationship
  Pair_Set pose_pairs;

  switch (i_mode)
  {
    case PAIR_MODE_EXHAUSTIVE:
      pose_pairs = exhaustivePairs(sfm_data.GetViews().size());
    break;
    case PAIR_MODE_CONTIGUOUS:
      pose_pairs = contiguousWithOverlap(vec_poses.size(), i_neighbor_count);
    break;
    case PAIR_MODE_NEIGHBORHOOD:
    {
      // List the poses priors
      std::vector<Vec3> vec_pose_centers;
      std::map<IndexT, IndexT> contiguous_to_pose_id;
      std::set<IndexT> used_pose_ids;
      for (const auto & view_it : sfm_data.GetViews())
      {
        const sfm::ViewPriors * prior = dynamic_cast<sfm::ViewPriors*>(view_it.second.get());
        if (prior != nullptr && prior->b_use_pose_center_ && used_pose_ids.count(prior->id_pose) == 0)
        {
          vec_pose_centers.push_back( prior->pose_center_ );
          contiguous_to_pose_id[contiguous_to_pose_id.size()] = prior->id_pose;
          used_pose_ids.insert(prior->id_pose);
        }
      }
      if (vec_pose_centers.empty())
      {
        std::cerr << "You are trying to use the gps_mode but your data does"
          << " not have any pose priors."
          << std::endl;
      }
      // Compute i_neighbor_count neighbor(s) for each pose
      matching::ArrayMatcherBruteForce<double> matcher;
      if (!matcher.Build(vec_pose_centers[0].data(), vec_pose_centers.size(), 3))
      {
        return EXIT_FAILURE;
      }
      size_t contiguous_pose_id = 0;
      for (const Vec3 pose_it : vec_pose_centers)
      {
        const double * query = pose_it.data();
        IndMatches vec_indices;
        std::vector<double> vec_distance;
        const int NN = i_neighbor_count + 1; // since itself will be found
        if (matcher.SearchNeighbours(query, 1, &vec_indices, &vec_distance, NN))
        {
          for (size_t i = 1; i < vec_indices.size(); ++i)
          {
            IndexT idxI = contiguous_to_pose_id.at(contiguous_pose_id);
            IndexT idxJ = contiguous_to_pose_id.at(vec_indices[i].j_);
            if (idxI > idxJ)
              std::swap(idxI, idxJ);
            pose_pairs.insert(Pair(idxI, idxJ));
          }
        }
        ++contiguous_pose_id;
      }
    }
    break;
    default:
      std::cerr << "Unknown pair mode." << std::endl;
      return EXIT_FAILURE;
  }


  // c. Convert the pose graph to a view graph
  Pair_Set view_pair;
  for (const auto & pose_pair : pose_pairs)
  {
    const IndexT poseA = pose_pair.first;
    const IndexT poseB = pose_pair.second;
    // get back the view related to those poses and create the pair (exhaustively)
    const auto range_a = pose_id_toViewId.equal_range(vec_poses[poseA]);
    for (auto view_id_a = range_a.first; view_id_a != range_a.second; view_id_a++)
    {
      const auto range_b = pose_id_toViewId.equal_range(vec_poses[poseB]);
      for (auto view_id_b = range_b.first; view_id_b != range_b.second; view_id_b++)
      {
        if (view_id_a != view_id_b)
        {
          view_pair.insert(
            Pair(std::min(view_id_a->second, view_id_b->second),
                 std::max(view_id_a->second, view_id_b->second)));
        }
      }
    }
  }

  if (view_pair.empty())
  {
    std::cout << "Warning: The computed pair list is empty...!" << std::endl;
  }

  // d. Export the view graph to a file and a SVG adjacency list

  AdjacencyMatrixToSVG(sfm_data.GetViews().size(), view_pair,
    stlplus::create_filespec(
      stlplus::folder_part(s_out_file),
      stlplus::filename_part(s_out_file), "svg"));

  if (savePairs(s_out_file, view_pair))
  {
    std::cout << "Exported " << view_pair.size() << " view pairs\n"
      <<"from a view graph that have " << pose_pairs.size()
      << " relative pose pairs." << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
