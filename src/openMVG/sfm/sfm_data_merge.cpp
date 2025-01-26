
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

bool badTrackRejector(double dPrecision, size_t count, SfM_Data& scene)
{
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
  const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(scene, dPrecision, 2);
    const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(scene, 2.0);

  return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
}

bool getVecs2Align(const openMVG::sfm::SfM_Data& first, const openMVG::sfm::SfM_Data& second, 
    std::vector<openMVG::Vec3> *parent_vecs,
    std::vector<openMVG::Vec3> *child_vecs,
    std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes)
{
  
  IndexT counter=0;
  std::vector<Pair> unused_cameras;

  for(auto p: sfm_filenames_indexes){
    std::pair<openMVG::IndexT, std::vector<openMVG::IndexT>> info = p.second;
    // lets make sure we have index in each SfM scene
    if(info.second.size()<2){continue;}
    
    IndexT p1 = info.second[0];
    IndexT p2 = info.second[1];

    //OPENMVG_LOG_INFO << p1 << "," << p2;

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
        parent_vecs -> push_back(first.poses.at(view1->id_pose).center());
        child_vecs -> push_back( second.poses.at(view2->id_pose).center() );
        counter++;
      }
    }
    catch(...){
      //std::cout << "Views not used in final reconstructions" << std::endl;
      unused_cameras.push_back(std::make_pair(p1,p2));
    }
  }

  if(parent_vecs->size() == child_vecs->size() && parent_vecs->size()>3){
    return true;
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
  //OPENMVG_LOG_INFO << "Non linear refinement" << std::endl;
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

bool mergeSfMScenes(openMVG::sfm::SfM_Data& sfm_data, openMVG::sfm::SfM_Data& child_sfm_data, 
   const double S, const openMVG::Mat3 R, const openMVG::Vec3 T,
   const std::map< std::string, std::pair<bool,std::vector<IndexT>> >& sfm_filenames_indexes
){

    std::set<Pair> common_ids;
    std::set<openMVG::IndexT> child_overlap_ids;

    for(auto pair: sfm_filenames_indexes){
        std::pair<openMVG::IndexT, std::vector<openMVG::IndexT>> info = pair.second;
        // lets make sure we have index in each SfM scene
        if(info.second.size()<2){continue;}
        
        IndexT p1 = info.second[0];
        IndexT p2 = info.second[1];

        common_ids.insert( Pair(p1,p2) );
        child_overlap_ids.insert(p2);
    }
    
    // references
    Views & parent_views = sfm_data.views;
    Poses & parent_poses = sfm_data.poses;

    Views & child_views = child_sfm_data.views;
    Poses & child_poses = child_sfm_data.poses;

    IndexT parent_views_size = parent_views.size();

    std::set<IndexT> remove_track_ids = std::set<IndexT>();

    std::set<Pair> new_view_pairings;
  
    for(auto & iterV : child_views){
        ViewPriors *prior = dynamic_cast<sfm::ViewPriors*>(iterV.second.get());
        IndexT id_view = prior->id_view;
        IndexT id_pose = prior->id_pose;

        // if the pose for the parent exists then we skip, if the pose doesn't we use 
        // the child scenes pose and update the parents
        if(child_overlap_ids.find(id_view) != child_overlap_ids.end()){

        for(auto p: common_ids){
            if(p.second!=id_view){continue;}
            const View * view1 = sfm_data.views.at(p.first).get();
            const View * view2 = child_sfm_data.views.at(p.second).get();

            if(!sfm_data.IsPoseAndIntrinsicDefined(view1) && child_sfm_data.IsPoseAndIntrinsicDefined(view2)){
                std::cout << "Pose reinstiated from child sfm scene " << std::endl;

                Pose3 pose = child_poses.at(view2->id_pose);
                Vec3 nloc = S * R * ( pose.center() ) + T; // update the camera position to the reference scene
                Pose3 npose = Pose3(pose.rotation(),nloc);

                parent_poses[view1->id_pose] = npose;

                break;
            }
        }

            continue;
        }

        // need to store the new view id to modify the observation ids
        new_view_pairings.insert(Pair(id_view,parent_views_size));

        // if there's a pose then modify the pose position
        if (child_sfm_data.IsPoseAndIntrinsicDefined(prior)){
            // image was not used in reconstruction
            // use the similarity transform and generate the new poses
            Pose3 pose = child_poses.at(id_pose);
            Vec3 nloc = S * R * ( pose.center() ) + T; // update the camera position to the reference scene
            Pose3 npose = Pose3(pose.rotation(),nloc);

            parent_poses[parent_views_size] = npose;

        }else{
            remove_track_ids.insert(id_view);
        }

        prior->id_view = parent_views_size;// new view id 
        prior->id_pose = parent_views_size;// new pose id

        parent_views[parent_views_size] = std::make_shared<ViewPriors>(*prior);

        parent_views_size++;
    }

    OPENMVG_LOG_INFO << "Views with no poses: " << remove_track_ids.size();

    Landmarks child_tracks = child_sfm_data.GetLandmarks();

    IndexT new_track_counter = sfm_data.structure.size();
    IndexT original_track_num = sfm_data.structure.size();
    // first update the landmarks in the child scene to the reference frame of the first
    for (auto& track: child_tracks)
    {
        IndexT track_id = track.first;
        Landmark landmark = track.second;
    
        Observations new_observations;
        for (auto& iterOb: landmark.obs)
        {
            IndexT id_view = iterOb.first;
            // need to update the view_id to the most up to date
            if(child_overlap_ids.find(id_view) != child_overlap_ids.end()){
                // check if the second view has any new information
                auto it = std::find_if(common_ids.begin(), common_ids.end(),
                [&](const Pair& val) -> bool {
                    return val.second == id_view;
                });
                id_view = it->first;
            }else{
                // else we want to use the first views pose if present
                auto it = std::find_if(new_view_pairings.begin(), new_view_pairings.end(),
                [&](const Pair& val) -> bool {
                    return val.first == id_view;
                });
                id_view = it->second;
            }

            try{
                const View * view = sfm_data.views.at(id_view).get();
                if(!sfm_data.IsPoseAndIntrinsicDefined(view)){
                    OPENMVG_LOG_WARNING << "Pose not defined for view " << id_view;
                    continue;
                }
            }catch(...){
                OPENMVG_LOG_WARNING << "View id does not exist " << id_view;
                continue;
            }
        

            // view id will have been updated by now, just have to insert it
            // iterOb.second.id_feat
            new_observations.insert({id_view, Observation(iterOb.second.x, UndefinedIndexT)});
            // new observations will have been added
        }

        Landmark new_landmark;
        new_landmark.X = S * R * ( landmark.X ) + T;
        new_landmark.obs = new_observations;

        sfm_data.structure[new_track_counter++] = std::move( new_landmark );
    }

    return true;
}

}//end of namespace sfm
}//end of namespace openmvg