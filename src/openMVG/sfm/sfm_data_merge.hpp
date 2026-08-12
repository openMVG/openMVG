#ifndef OPENMVG_SFM_SFM_MERGE_HPP
#define OPENMVG_SFM_SFM_MERGE_HPP

#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#include "openMVG/system/logger.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "ceres/problem.h"
#include "ceres/solver.h"

namespace openMVG {
namespace sfm {

// filename -> (is_new_to_first, [view id in first (if any), view id in second (if any)])
using OverlapMap = std::unordered_map<std::string, std::pair<bool, std::vector<IndexT>>>;

std::set<IndexT> GetValidIntrinsicsIds(const SfM_Data& sfm_data1);

std::set<IndexT> getCommonCameraIds(const SfM_Data& sfm_data_1, const SfM_Data& sfm_data_2);

std::shared_ptr<cameras::IntrinsicBase> findBestIntrinsic(const SfM_Data& sfm_data_A, const SfM_Data& sfm_data_B, openMVG::IndexT cam_id);

/*
*   Merge every camera intrinsic from second_sfm_data into sfm_data.
*   Intrinsic ids that exist in both scenes are assumed to describe the same physical
*   camera: the better of the two (by reprojection RMSE, via findBestIntrinsic) is kept.
*   Intrinsic ids that only exist in second_sfm_data are copied in under a freshly
*   allocated id.
*   sfm_data.intrinsics is updated in place.
*   @return a map from second_sfm_data intrinsic id -> merged (sfm_data) intrinsic id,
*           to be used to remap View::id_intrinsic when copying views across.
*/
std::unordered_map<IndexT, IndexT> mergeIntrinsics(SfM_Data& sfm_data, const SfM_Data& second_sfm_data);

std::vector<std::string> newUniqueImages(const OverlapMap& sfm_filenames_indexes);

// count the number of images that are duplicates
openMVG::IndexT getNumberOverlapping(const OverlapMap& sfm_filenames_indexes);

// print information of the overlap
std::string printOverlapInformation(const OverlapMap& sfm_filenames_indexes);

/*
*   This function gets the overlapping pairs of image names between 2 SfM scenes
*   it returns bool if the scenes are Overlapping and the overlapping_pairs std::map
*   contains a tuple with first element being the image name, index in scene 1 , index in scene 2
*/
bool getOverlappingImages(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second,
    OverlapMap& sfm_filenames_indexes);

bool badTrackRejector(double dPrecision, size_t count, SfM_Data& scene);

bool getVecs2Align(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second,
    std::vector<openMVG::Vec3> *first_vecs,
    std::vector<openMVG::Vec3> *second_vecs,
    OverlapMap& sfm_filenames_indexes);

bool computeSimilarity(
  const std::vector<openMVG::Vec3> & vec_camPosGT,
  const std::vector<openMVG::Vec3> & vec_camPosComputed,
  std::vector<openMVG::Vec3> & vec_camPosComputed_T,
  double *Sout, openMVG::Mat3 * Rout, openMVG::Vec3 * tout);

bool mergeSfMScenes(openMVG::sfm::SfM_Data& sfm_data, openMVG::sfm::SfM_Data& second_sfm_data,
   const double S, const openMVG::Mat3 R, const openMVG::Vec3 T,
   const OverlapMap& sfm_filenames_indexes,
   const std::unordered_map<IndexT, IndexT>& intrinsic_id_remap,
   double landmark_fusion_pixel_tolerance = 4.0);

}//end of namespace sfm
}//end of namespace openmvg

#endif // OPENMVG_SFM_SFM_MERGE_HPP
