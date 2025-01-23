#ifndef OPENMVG_SFM_SFM_MERGE_HPP
#define OPENMVG_SFM_SFM_MERGE_HPP

#include <vector>
#include <algorithm>
#include <set>
#include <iostream>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"

#include "openMVG/system/logger.hpp"

namespace openMVG {
namespace sfm {

std::set<IndexT> GetValidIntrinsicsIds(const SfM_Data& sfm_data1);

std::set<IndexT> getCommonCameraIds(const SfM_Data& sfm_data_1, const SfM_Data& sfm_data_2);

std::shared_ptr<cameras::IntrinsicBase> findBestIntrinsic(const SfM_Data& sfm_data_A, const SfM_Data& sfm_data_B, openMVG::IndexT cam_id);

std::vector<std::string> newUniqueImages(const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes);

// count the number of images that are duplicates
openMVG::IndexT getNumberOverlapping(const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes);

// print information of the overlap 
std::string printOverlapInformation(const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes);

/*
*   This function gets the overlapping pairs of image names between 2 SfM scenes
*   it returns bool if the scenes are Overlapping and the overlapping_pairs std::map
*   contains a tuple with first element being the image name, index in scene 1 , index in scene 2
*/
bool getOverlappingImages(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second, 
    std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes);

bool getVecs2Align(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second, 
    std::vector<openMVG::Vec3>& parent_vecs,
    std::vector<openMVG::Vec3>& child_vecs,
    std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes);

bool computeSimilarity(
  const std::vector<openMVG::Vec3> & vec_camPosGT,
  const std::vector<openMVG::Vec3> & vec_camPosComputed,
  std::vector<openMVG::Vec3> & vec_camPosComputed_T,
  double *Sout, openMVG::Mat3 * Rout, openMVG::Vec3 * tout);

bool mergeSfMScenes(openMVG::sfm::SfM_Data& sfm_data, openMVG::sfm::SfM_Data& child_sfm_data, 
   const double S, const openMVG::Mat3 R, const openMVG::Vec3 T,
   const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes);

}//end of namespace sfm
}//end of namespace openmvg

#endif // OPENMVG_SFM_SFM_MERGE_HPP