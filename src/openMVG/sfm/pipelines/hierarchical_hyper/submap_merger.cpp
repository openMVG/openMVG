// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "submap_merger.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/sfm_data_BA_fixed_points.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_utilities.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "ceres/types.h"

namespace openMVG{
namespace sfm{

SubmapMerger::SubmapMerger(const HsfmSubmaps & submaps, const cameras::Intrinsic_Parameter_Type & intrinsic_refinement_options)
  : submaps_(submaps), scene_aligner_(new SceneAligner), intrinsic_refinement_options_(intrinsic_refinement_options)
{}

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
      const HsfmSubmap & smapA = submaps_.at(*smap_id);
      const openMVG::IndexT parent_id = smapA.parent_id;

      const IndexT sister_smap_id = getSiblingSubmapId(submaps_, *smap_id);

      // check that sister submap is in the "to merge" list. if not, continue to next submap.
      if (submaps_to_merge.find(sister_smap_id) == submaps_to_merge.end())
      {
        ++smap_id;
        continue;
      }
      std::cout << "Merging submap " << *smap_id << " with submap " << sister_smap_id << std::endl;

      // Merge the two submaps together
      MergeSubmapPair(parent_id, debug_out_dir);

      // remove indices from the to-merge list !
      submaps_to_merge.erase(sister_smap_id);
      smap_id = submaps_to_merge.erase(smap_id);

      // add parent submap to the to-merge list !
      submaps_to_merge.insert(parent_id);
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

  // --- OPTIMIZATION OF THE BASE NODE POSES + SEPARATOR LANDMARKS POSITIONS ---
  //
  // We have to find the best similarity transform between the first and the second
  // submap. This transformation is a classic combination of translation (3 parameters),
  // rotation (3 parameters, ceres uses angle-axis representation), and scaling (1 parameter).
  // Using a 3d Bundle adjustment step, we refine these 7 parameters and find the relative
  // pose+scaling of the second submap with respect to the first one.
  // For this optimization, we use the landmarks that are common to both submaps (a subset of
  // the separators landmarks), and minimize the error of the distance
  // between the separator landmarks and the base node positions.

  Bundle_Adjustment_Ceres::BA_Ceres_options options;
  options.linear_solver_type_ = ceres::DENSE_SCHUR;

  if (!MergeScenesUsingCommonTracks(parent_sfm_data, first_sfm_data, second_sfm_data, parent_submap.separator, scene_aligner_.get()))
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

  std::cout << "\n Optimize newly merged parent submap ! \n";
  // --- ADDITIONAL OPTIMIZATION OF THE NEWLY FILLED PARENT SUBMAP ---
  // This is like a normal bundle adjustment step of the whole parent
  // submap where we keep the position of the separator landmarks as FIXED. (as advised in master thesis [TODO ref])

  // note : parameters copied from sequential sfm
  const double requiredPixelResidualError = 4.0;
  const size_t outlierNumberThreshold = 50;

  do
  {
    BundleAdjustment_FixedSeparators(parent_id);
  }
  while (badTrackRejector(requiredPixelResidualError, outlierNumberThreshold, parent_sfm_data));
  eraseUnstablePosesAndObservations(parent_sfm_data);

  return true;
}

/**
 * @brief Discard tracks with too large residual error
 *
 * Remove observation/tracks that have:
 *  - too large residual error
 *  - too small angular value
 *
 * @note copied over from sequential sfm !
 * @return True if more than 'count' outliers have been removed.
 */
bool SubmapMerger::badTrackRejector(double dPrecision, size_t count, SfM_Data& scene)
{
  const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(scene, dPrecision, 2);
  const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(scene, 2.0);

  return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
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
  const std::unique_ptr<Bundle_Adjustment_Fixed_Points> bundle_adjustment_ptr =
      std::unique_ptr<Bundle_Adjustment_Fixed_Points>(new Bundle_Adjustment_Fixed_Points(options));

  const Optimize_Options ba_refine_options
    ( intrinsic_refinement_options_,
      Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
      Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
      Control_Point_Parameter(),
      false // no motion prior for now... TODO add this ?
    );
  return bundle_adjustment_ptr->Adjust(sfm_data, submap.separator, ba_refine_options);
}

}
}
