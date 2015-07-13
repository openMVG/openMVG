
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "openMVG/stl/stl.hpp"

#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace sfm {

using namespace openMVG::geometry;
using namespace openMVG::cameras;
using namespace openMVG::image;


/// Find the color of the SfM_Data Landmarks/structure
void ColorizeTracks( SfM_Data & sfm_data )
{
  // Colorize each track
  //  Start with the most representative image
  //    and iterate to provide a color to each 3D point

  {
    std::vector<Vec3> vec_3dPoints;
    std::vector<Vec3> vec_tracksColor;

    C_Progress_display my_progress_bar(sfm_data.GetLandmarks().size(),
                                       std::cout,
                                       "\nCompute scene structure color\n");

    vec_3dPoints.resize(sfm_data.GetLandmarks().size());

    //Build a list of contiguous index for the trackIds
    std::map<IndexT, IndexT> trackIds_to_contiguousIndexes;
    IndexT cpt = 0;
    for (Landmarks::const_iterator it = sfm_data.GetLandmarks().begin();
      it != sfm_data.GetLandmarks().end(); ++it, ++cpt)
    {
      trackIds_to_contiguousIndexes[it->first] = cpt;
      vec_3dPoints[cpt] = it->second.X;
    }

    // The track list that will be colored (point removed during the process)
    std::set<IndexT> remainingTrackToColor;
    std::transform(sfm_data.GetLandmarks().begin(), sfm_data.GetLandmarks().end(),
      std::inserter(remainingTrackToColor, remainingTrackToColor.begin()),
      stl::RetrieveKey());

    while( !remainingTrackToColor.empty() )
    {
      // Find the most representative image (for the remaining 3D points)
      //  a. Count the number of observation per view for each 3Dpoint Index
      //  b. Sort to find the most representative view index

      std::map<IndexT, IndexT> map_IndexCardinal; // ViewId, Cardinal
      for (std::set<IndexT>::const_iterator
        iterT = remainingTrackToColor.begin();
        iterT != remainingTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = *iterT;
        const Observations & obs = sfm_data.GetLandmarks().at(trackId).obs;
        for( Observations::const_iterator iterObs = obs.begin();
          iterObs != obs.end(); ++iterObs)
        {
          const size_t viewId = iterObs->first;
          if (map_IndexCardinal.find(viewId) == map_IndexCardinal.end())
            map_IndexCardinal[viewId] = 1;
          else
            ++map_IndexCardinal[viewId];
        }
      }

      // Find the View index that is the most represented
      std::vector<IndexT> vec_cardinal;
      std::transform(map_IndexCardinal.begin(),
        map_IndexCardinal.end(),
        std::back_inserter(vec_cardinal),
        stl::RetrieveValue());
      using namespace stl::indexed_sort;
      std::vector< sort_index_packet_descend< IndexT, IndexT> > packet_vec(vec_cardinal.size());
      sort_index_helper(packet_vec, &vec_cardinal[0], 1);

      // First image index with the most of occurence
      std::map<IndexT, IndexT>::const_iterator iterTT = map_IndexCardinal.begin();
      std::advance(iterTT, packet_vec[0].index);
      const size_t view_index = iterTT->first;
      const View * view = sfm_data.GetViews().at(view_index).get();
      const std::string sView_filename = stlplus::create_filespec(sfm_data.s_root_path,
        view->s_Img_path);
      Image<RGBColor> image;
      ReadImage(sView_filename.c_str(), &image);

      // Iterate through the remaining track to color
      // - look if the current view is present to color the track
      std::set<IndexT> set_toRemove;
      for (std::set<IndexT>::const_iterator
        iterT = remainingTrackToColor.begin();
        iterT != remainingTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = *iterT;
        const Observations & obs = sfm_data.GetLandmarks().at(trackId).obs;
        Observations::const_iterator it = obs.find(view_index);

        if (it != obs.end())
        {
          // Color the track
          const Vec2 & pt = it->second.x;
          sfm_data.structure.at(trackId).rgb = image(pt.y(), pt.x());
          set_toRemove.insert(trackId);
          ++my_progress_bar;
        }
      }
      // Remove colored track
      for (std::set<IndexT>::const_iterator iter = set_toRemove.begin();
        iter != set_toRemove.end(); ++iter)
      {
        remainingTrackToColor.erase(*iter);
      }
    }
  }
}

} // namespace sfm
} // namespace openMVG