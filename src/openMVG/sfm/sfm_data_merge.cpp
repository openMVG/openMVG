
// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_merge.hpp"

namespace openMVG {
namespace sfm {

/**
* @brief Get the best intrinsics for this particular intrinsic id
* @return the lowest RMSE camera intrinsic shared_ptr object
* @note return is based on the best RMSE from that partiular intrinsic id
*/
std::shared_ptr<cameras::IntrinsicBase> findBestIntrinsic(const SfM_Data& sfm_data_A, const SfM_Data& sfm_data_B, openMVG::IndexT cam_id)
{
  // both alternatives for the intrinsic
  const auto intrinsicA = sfm_data_A.intrinsics.at(cam_id);
  const auto intrinsicB = sfm_data_B.intrinsics.at(cam_id);

  // quick escape in case intrinsics are equal
  if (intrinsicA == intrinsicB // check pointer equality
      || (intrinsicA->getType() == intrinsicB->getType() // check intrinsics equality
         && intrinsicA->hashValue() == intrinsicB->hashValue())){
    return std::shared_ptr<cameras::IntrinsicBase>(intrinsicA->clone());
  }
    

  double RMSE_A(0.0), RMSE_B(0.0);
  int n_totalResiduals(0);
  for (const auto & sfm_data : {sfm_data_A, sfm_data_B})
  {
    for (const auto & landmark : sfm_data.GetLandmarks())
    {
      const Observations & observations = landmark.second.obs;
      for (const auto & obs: observations)
      {
        // we have to do the following check because observations of common landmarks are not pruned out when
        // clustering a submap in two...which makes it simpler to merge back together
        const auto & view = sfm_data.GetViews().find(obs.first);
        if (view == sfm_data.GetViews().end()){
          continue;
        }

        const IndexT & intrinsic_id = sfm_data.GetViews().at(obs.first)->id_intrinsic;
        if (intrinsic_id != cam_id){
          continue;
        }

        const IndexT & pose_id = sfm_data.GetViews().at(obs.first)->id_pose;
        const auto & pose = sfm_data.GetPoses().at(pose_id);
        const Vec3 X = pose(landmark.second.X);
        const Vec2 residualA = intrinsicA->residual(X, obs.second.x);
        const Vec2 residualB = intrinsicB->residual(X, obs.second.x);
        RMSE_A += residualA(0) * residualA(0);
        RMSE_A += residualA(1) * residualA(1);
        RMSE_B += residualB(0) * residualB(0);
        RMSE_B += residualB(1) * residualB(1);
        ++n_totalResiduals;
      }
    }
  }

  RMSE_A = std::sqrt((RMSE_A)/(n_totalResiduals));
  RMSE_B = std::sqrt((RMSE_B)/(n_totalResiduals));
  std::cout << "RMSE_A : " << RMSE_A << " RMSE_B : "  << RMSE_B << std::endl;

  if (RMSE_A < RMSE_B){
    return std::shared_ptr<cameras::IntrinsicBase>(intrinsicA->clone());
  }
  else{
    return std::shared_ptr<cameras::IntrinsicBase>(intrinsicB->clone());
  }
}

/**
* @brief Get a list of valid intrinsic Ids from sfm file
* @return returns a list of valid intrinsics in the sfm_data file
* @note returns the indexes of the valid intrinsics
*/
std::set<IndexT> GetValidIntrinsicsIds(const SfM_Data& sfm_data1){
  std::set<IndexT> valid_views;
  for (const auto & view : sfm_data1.GetViews())
  {
    const View * v = view.second.get();
    if (sfm_data1.GetIntrinsics().find(v->id_intrinsic) != sfm_data1.GetIntrinsics().end()){
      valid_views.insert(v->id_view);
    }
  }
  return valid_views;
}

/**
* @brief Get a list of valid intrinsic Ids common in both sfm scenes
* @return returns a list of valid intrinsics in both sfm_data files
* @note if not common in both, it'll be left out
*/
std::set<IndexT> getCommonCameraIds(const SfM_Data& sfm_data_1, const SfM_Data& sfm_data_2)
{
  const std::set<IndexT> cam_ids_1 = GetValidIntrinsicsIds(sfm_data_1);
  const std::set<IndexT> cam_ids_2 = GetValidIntrinsicsIds(sfm_data_2);

  std::set<IndexT> common_cam_ids;
  std::set_intersection(cam_ids_1.cbegin(), cam_ids_1.cend(),
                        cam_ids_2.cbegin(), cam_ids_2.cend(),
                        std::inserter(common_cam_ids, common_cam_ids.begin()));
  return common_cam_ids;
}

std::vector<std::string> newUniqueImages(const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes){
    std::vector<std::string> new_filenames;
    for(auto p: sfm_filenames_indexes){
        std::pair<openMVG::IndexT, std::vector<openMVG::IndexT>> info = p.second;
        if(info.first){new_filenames.push_back(p.first);}
    }
    return new_filenames;
}

/**
* @brief See how many images are overlapping
* @return returns a single IndexT value, 0 for no Overlapping, >0 for overlap
* @note returns the indexes of the valid intrinsics
*/
openMVG::IndexT getNumberOverlapping(const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes){
    openMVG::IndexT counter = 0;
    for(auto p: sfm_filenames_indexes){
        std::pair<openMVG::IndexT, std::vector<openMVG::IndexT>> info = p.second;
        if(info.second.size() > 1){counter++;}
    }
    return counter;
}

/**
* @brief See how many images are overlapping
* @return returns a single IndexT value, 0 for no Overlapping, >0 for overlap
* @note returns the indexes of the valid intrinsics
*/
std::string printOverlapInformation(const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes){
    
    std::string result = "Name, Size, New Image, Indexes\n";
    
    for(auto p: sfm_filenames_indexes){
        std::pair<openMVG::IndexT, std::vector<openMVG::IndexT>> info = p.second;

        result += p.first+","+std::to_string(info.second.size())+",";
        std::string indexes = " ";
        for(auto a: info.second){
            indexes += std::to_string(a) + " ";
        }
        result +=  indexes + "\n";
    }
    return result;
}

/**
* @brief Get a map of all filenames in 2 sfm scenes, each filename is key, then second is a vector of indexes in first and second
* @return return a map of std::string and std::vector, filename key, vector contains indexes in 1 and 2
* @note returns the indexes of the valid intrinsics
*/
bool getOverlappingImages(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second, 
    std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes){
  // the index and strings will all be different is the thinking here

  IndexT duplicates = 0;

  for (auto iterViews = first.views.cbegin();
        iterViews != first.views.cend();++iterViews){
      const openMVG::sfm::View * view = iterViews->second.get();
      std::string filename = view->s_Img_path;
      IndexT index = view->id_view;
      sfm_filenames_indexes.insert( std::make_pair(filename, std::make_pair(false, std::vector<IndexT>{index})) );
  }

  for (auto iterViews = second.views.cbegin();
        iterViews != second.views.cend();++iterViews){
      const openMVG::sfm::View * view = iterViews->second.get();
      std::string filename = view->s_Img_path;
      IndexT index = view->id_view;
      
      if(sfm_filenames_indexes.count(filename) > 0){
        sfm_filenames_indexes[filename].second.push_back(index);
        duplicates++;
      }else{
        sfm_filenames_indexes.insert( std::make_pair(filename, std::make_pair(false, std::vector<IndexT>{index})) );
      }
  }
  if(duplicates>0){return true;}
  return false;
}

bool getVecs2Align(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second, 
    openMVG::IndexT& overlap_amount, std::vector<openMVG::Vec3>& parent_vecs,
    std::vector<openMVG::Vec3>& child_vecs,
    std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes)
{
  
  IndexT counter=0;
  std::vector<Pair> unused_cameras;

  for(auto p: sfm_filenames_indexes){
    std::pair<openMVG::IndexT, std::vector<openMVG::IndexT>> info = p.second;
    if(info.second.size()<1){continue;}

    IndexT p1 = info.second[0];
    IndexT p2 = info.second[1];
    try{
      const View * view1 = first.views.at(p1).get();
      const View * view2 = second.views.at(p2).get();

      if(!first.IsPoseAndIntrinsicDefined(view1)){
        OPENMVG_LOG_INFO << view1->id_view << " has no pose associated in parent";
        unused_cameras.push_back(std::make_pair(p1,p2));
      }
      else if(!second.IsPoseAndIntrinsicDefined(view2)){
        OPENMVG_LOG_INFO << view2->id_view << " has no pose associated in child";
        unused_cameras.push_back(std::make_pair(p1,p2));
      }
      else{
        parent_vecs[counter] = first.poses.at(view1->id_pose).center();
        child_vecs[counter] = second.poses.at(view2->id_pose).center();
        counter++;
      }
    }
    catch(...){
      //std::cout << "Views not used in final reconstructions" << std::endl;
      unused_cameras.push_back(std::make_pair(p1,p2));
    }
  }

  return false;
}

/// Compute a 7DOF rigid transform between the two camera trajectories
bool computeSimilarity(
  const std::vector<openMVG::Vec3> & vec_camPosGT,
  const std::vector<openMVG::Vec3> & vec_camPosComputed,
  std::vector<openMVG::Vec3> & vec_camPosComputed_T,
  double *Sout, openMVG::Mat3 * Rout, openMVG::Vec3 * tout)
{
  if (vec_camPosGT.size() != vec_camPosComputed.size()) {
    OPENMVG_LOG_ERROR << "Cannot perform registration, vector sizes are different";
    return false;
  }

  // Move input point in appropriate container
  openMVG::Mat x1(3, vec_camPosGT.size());
  openMVG::Mat x2(3, vec_camPosGT.size());
  for (size_t i = 0; i  < vec_camPosGT.size(); ++i) {
    x1.col(i) = vec_camPosComputed[i];
    x2.col(i) = vec_camPosGT[i];
  }
  // Compute rigid transformation p'i = S R pi + t

  double S;
  Vec3 t;
  Mat3 R;
  openMVG::geometry::FindRTS(x1, x2, &S, &t, &R);
  OPENMVG_LOG_INFO << "Non linear refinement" << std::endl;
  openMVG::geometry::Refine_RTS(x1,x2,&S,&t,&R);

  vec_camPosComputed_T.resize(vec_camPosGT.size());
  for (size_t i = 0; i  < vec_camPosGT.size(); ++i)
  {
    const Vec3 newPos = S * R * ( vec_camPosComputed[i]) + t;
    vec_camPosComputed_T[i] = newPos;
  }

  *Sout = S;
  *Rout = R;
  *tout = t;
  return true;
}

}//end of namespace sfm
}//end of namespace openmvg