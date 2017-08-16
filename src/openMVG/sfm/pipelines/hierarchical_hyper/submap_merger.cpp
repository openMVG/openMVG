#include "submap_merger.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_BA_fixed_separators.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_utilities.hpp"

#include "ceres/types.h"

namespace openMVG{
namespace sfm{

std::vector<size_t> findCommonReconstructedSeparatorTracks(HsfmSubmap& parent_submap, const openMVG::sfm::SfM_Data& first_sfm_data, const openMVG::sfm::SfM_Data& second_sfm_data)
{
  std::vector<size_t> common_track_ids;
  const Landmarks & first_landmarks = first_sfm_data.GetLandmarks();
  const Landmarks & second_landmarks = second_sfm_data.GetLandmarks();
  std::cout << first_landmarks.size() << " " << second_landmarks.size() << std::endl;
  for (const auto & track_id : parent_submap.separator)
  {
    if (first_landmarks.find(track_id) == first_landmarks.end()
        || second_landmarks.find(track_id) == second_landmarks.end())
      continue;

    // if the landmark is reconstructed in both submaps -> add it to the common tracks
    common_track_ids.push_back(track_id);
  }
  std::cout << common_track_ids.size() << " tracks in common" << std::endl;
  return common_track_ids;
}

SubmapMerger::SubmapMerger(const HsfmSubmaps & submaps)
  : submaps_(submaps)
{
  scene_aligner_ = std::unique_ptr<SceneAligner>(new SceneAligner);
}

bool SubmapMerger::Merge(const std::string & debug_out_dir)
{
  // - loop over leaf submaps
  // - find corresponding sister submap
  // - merge with sister submap

  // store ids of submaps to merge together, initialize this list with the leaves submaps.
  std::set<IndexT> submaps_to_merge; // "to-merge" list
  for (const auto & smap : submaps_)
  {
    if (smap.second.is_parent)
      continue;
    submaps_to_merge.insert(smap.first);
    std::cout << smap.first << std::endl;
  }

  // as long as there are submaps left to merge, keep merging
  // NOTE : by convention, submap 0 is the root submap. We thus stop this loop when
  // submap 0 is the only submap left.
  while (!submaps_to_merge.empty() && submaps_to_merge.count(0) == 0)
  {
    // Loop over submaps left to merge
    for (std::set<IndexT>::iterator smap_id = submaps_to_merge.begin(); smap_id != submaps_to_merge.end();)
    {
      HsfmSubmap & smap = submaps_.at(*smap_id);

      const IndexT sister_smap_id = getSiblingSubmapId(submaps_, *smap_id);

      // check that sister submap is in the "to merge" list. if not, continue to next submap.
      if (submaps_to_merge.find(sister_smap_id) == submaps_to_merge.end())
      {
        smap_id++;
        continue;
      }

      std::cout << "Merging submap " << *smap_id << " with submap " << sister_smap_id << std::endl;

      // Merge the two submaps together
      MergeSubmapPair(smap.parent_id, debug_out_dir);

      // remove indices from the to-merge list !
      submaps_to_merge.erase(sister_smap_id);
      smap_id = submaps_to_merge.erase(smap_id);

      // add parent submap to the to-merge list !
      submaps_to_merge.insert(smap.parent_id);
    }
  }
  return true;
}

bool SubmapMerger::MergeSubmapPair(const IndexT parent_id, const std::string &debug_out_dir)
{
  HsfmSubmap & parent_submap = submaps_.at(parent_id);
  SfM_Data & parent_sfm_data = parent_submap.sfm_data;
  const HsfmSubmap & first_submap = submaps_.at(submaps_.at(parent_id).children_submaps.first);
  const HsfmSubmap & second_submap = submaps_.at(submaps_.at(parent_id).children_submaps.second);
  const openMVG::sfm::SfM_Data & first_sfm_data = first_submap.sfm_data;
  const openMVG::sfm::SfM_Data & second_sfm_data = second_submap.sfm_data;

  // find all common reconstructed tracks in both submaps
  // (note that this is a subset of the separator tracks)
  std::vector<size_t> common_track_ids = findCommonReconstructedSeparatorTracks(parent_submap, first_sfm_data, second_sfm_data);

  // --- OPTIMIZATION OF THE BASE NODE POSES + SEPARATOR LANDMARKS POSITIONS ---
  //
  // We have to find the best transformation between the first and the second
  // submap. This transformation is a classic combination of translation (3 parameters),
  // rotation (3 parameters, ceres uses angle-axis representation), and scaling (1 parameter).
  // Using a 3d Bundle adjustment step, we refine these 7 parameters and find the relative
  // pose+scaling of the second submap with respect to the first one.
  // For this optimization, we use the landmarks that are common to both submaps (a subset of
  // the separators landmarks, contained in common_track_ids), and minimize the error of the distance
  // between the separator landmarks and the base node positions.

  Bundle_Adjustment_Ceres::BA_Ceres_options options;
  options.linear_solver_type_ = ceres::DENSE_SCHUR;

  if (!MergeScenesUsingCommonTracks(parent_sfm_data, first_sfm_data, second_sfm_data, common_track_ids, scene_aligner_.get()))
  {
    std::cerr << "Submaps merging failed !" << std::endl;
    return false;
  }

  // DEBUG - export submap before optimization
  if (debug_out_dir != "")
  {
    ExportSubmapData(submaps_, parent_id, stlplus::create_filespec(debug_out_dir, "sfm_data_" + std::to_string(parent_id) + "_beforeopt", "ply"));
  }

  // TODO : (could be done), add eventual camera poses that were not reconstructed in
  // the submaps due to losing some information while clustering ? (i.e. a few more resections
  // in the parent submap with all the data from the two children submaps).

  // --- ADDITIONAL OPTIMIZATION OF THE NEWLY FILLED PARENT SUBMAP ---
  // This is like a normal bundle adjustment step of the whole parent
  // submap where we keep the position of the separator landmarks as FIXED. (as advised in master thesis [TODO ref])

  std::cout << "\n Optimize newly merged parent submap ! \n";
  BundleAdjustment_FixedSeparators(parent_id);

  return true;
}

bool SubmapMerger::BundleAdjustment_FixedSeparators(const IndexT smap_id)
{
  HsfmSubmap & submap = submaps_.at(smap_id);
  SfM_Data & sfm_data = submap.sfm_data;

  Bundle_Adjustment_Ceres::BA_Ceres_options options;
  if ( sfm_data.GetPoses().size() > 100 &&
      (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
      )
  // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
  {
    options.preconditioner_type_ = ceres::JACOBI;
    options.linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
  }
  const std::unique_ptr<Bundle_Adjustment_Fixed_Separators> bundle_adjustment_ptr =
      std::unique_ptr<Bundle_Adjustment_Fixed_Separators>(new Bundle_Adjustment_Fixed_Separators(options));
  return bundle_adjustment_ptr->Adjust(sfm_data, submap.separator);
}

}
}
