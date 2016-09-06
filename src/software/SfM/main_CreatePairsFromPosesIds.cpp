// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/system/timer.hpp"

using namespace openMVG;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/progress/progress.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <map>

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Given a linear pose graph and an given overlap value:\n"
    << " - export corresponding view pairs file that can be used in \n"
    << "     main_ComputeMatches.\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string s_SfM_Data_filename;
  std::string s_out_file;
  int i_overlapping = 5;

  cmd.add( make_option('i', s_SfM_Data_filename, "input_file") );
  cmd.add( make_option('o', s_out_file, "out_file") );
  cmd.add( make_option('p', i_overlapping, "pose_overlapping") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-o|--out_file] the output pairlist file\n"
    << "optional:\n"
    << "[-p|--pose_overlapping] Define the pose neighborhood for matching\n"
    << "\t(default: " << i_overlapping << ")\n"
    << "\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, s_SfM_Data_filename, ESfM_Data(VIEWS|INTRINSICS)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << s_SfM_Data_filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

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
  // a. List the view pose
  // b. Establish a pose graph 'linear pose relation' with a bounded link coverage
  // c. Convert the pose graph edges to a view graph
  //---------------------------------------
  std::multimap<IndexT, IndexT> pose_id_toViewId;
  std::set<IndexT> set_poses;

  // a. Get nodes of the pose graph
  for (const auto & viewIter : sfm_data.GetViews())
  {
    const View * v = viewIter.second.get();
    assert (viewIter.first == v->id_view);
    pose_id_toViewId.insert( std::make_pair(v->id_pose, v->id_view) );
    set_poses.insert(v->id_pose);
  }
  const std::vector<IndexT> vec_poses(set_poses.begin(), set_poses.end());

  // b. Create the 'linear' pose graph pair relationship
  const Pair_Set pose_pairs = contiguousWithOverlap(set_poses.size(), i_overlapping);

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

  if (savePairs(s_out_file, view_pair))
  {
    std::cout << "Exported " << view_pair.size() << " view pairs\n"
      <<"from a view graph that have " << pose_pairs.size()
      << " relative pose pairs." << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}

