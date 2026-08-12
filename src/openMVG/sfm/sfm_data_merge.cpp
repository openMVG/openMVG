// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_merge.hpp"

#include <numeric>

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
  // iterate by pointer: {sfm_data_A, sfm_data_B} would otherwise build a temporary
  // std::initializer_list<SfM_Data>, deep-copying both entire scenes (views, poses,
  // structure, intrinsics) just to read them.
  for (const SfM_Data * sfm_data : {&sfm_data_A, &sfm_data_B})
  {
    for (const auto & landmark : sfm_data->GetLandmarks())
    {
      const Observations & observations = landmark.second.obs;
      for (const auto & obs: observations)
      {
        // we have to do the following check because observations of common landmarks are not pruned out when
        // clustering a submap in two...which makes it simpler to merge back together
        const auto & view = sfm_data->GetViews().find(obs.first);
        if (view == sfm_data->GetViews().end()){
          continue;
        }

        const IndexT & intrinsic_id = sfm_data->GetViews().at(obs.first)->id_intrinsic;
        if (intrinsic_id != cam_id){
          continue;
        }

        const IndexT & pose_id = sfm_data->GetViews().at(obs.first)->id_pose;
        const auto & pose = sfm_data->GetPoses().at(pose_id);
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
  OPENMVG_LOG_INFO << "RMSE_A : " << RMSE_A << " RMSE_B : "  << RMSE_B ;

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

/**
* @brief Merge every camera intrinsic from second_sfm_data into sfm_data
* @note intrinsic ids present in both scenes are assumed to be the same physical camera:
*       the better of the two (by findBestIntrinsic RMSE) is kept under that id.
*       intrinsic ids only present in second_sfm_data are copied in under a fresh id.
* @return map from second_sfm_data intrinsic id -> merged (sfm_data) intrinsic id
*/
std::unordered_map<IndexT, IndexT> mergeIntrinsics(SfM_Data& sfm_data, const SfM_Data& second_sfm_data)
{
  std::unordered_map<IndexT, IndexT> intrinsic_id_remap;

  IndexT next_intrinsic_id = 0;
  for (const auto & intrinsic : sfm_data.GetIntrinsics())
  {
    next_intrinsic_id = std::max(next_intrinsic_id, intrinsic.first + 1);
  }

  for (const auto & second_intrinsic : second_sfm_data.GetIntrinsics())
  {
    const IndexT second_id = second_intrinsic.first;

    if (sfm_data.intrinsics.find(second_id) != sfm_data.intrinsics.end())
    {
      // Both scenes define an intrinsic under this id: treat them as the same physical
      // camera and keep whichever one better explains the observed residuals.
      sfm_data.intrinsics[second_id] = findBestIntrinsic(sfm_data, second_sfm_data, second_id);
      intrinsic_id_remap[second_id] = second_id;
    }
    else
    {
      // A camera only present in the second scene: bring it in under a fresh id.
      const IndexT new_id = next_intrinsic_id++;
      sfm_data.intrinsics[new_id] = std::shared_ptr<cameras::IntrinsicBase>(second_intrinsic.second->clone());
      intrinsic_id_remap[second_id] = new_id;
    }
  }

  return intrinsic_id_remap;
}

std::vector<std::string> newUniqueImages(const OverlapMap& sfm_filenames_indexes){
    std::vector<std::string> new_filenames;
    for(const auto & p: sfm_filenames_indexes){
        const std::pair<bool,std::vector<IndexT>> & info = p.second;
        if(info.first){new_filenames.push_back(p.first);}
    }
    return new_filenames;
}

/**
* @brief See how many images are overlapping
* @return returns a single IndexT value, 0 for no Overlapping, >0 for overlap
* @note returns the indexes of the valid intrinsics
*/
openMVG::IndexT getNumberOverlapping(const OverlapMap& sfm_filenames_indexes){
    openMVG::IndexT counter = 0;
    for(const auto & p: sfm_filenames_indexes){
        const std::pair<bool,std::vector<IndexT>> & info = p.second;
        if(info.second.size() > 1){counter++;}
    }
    return counter;
}

/**
* @brief See how many images are overlapping
* @return returns a single IndexT value, 0 for no Overlapping, >0 for overlap
* @note returns the indexes of the valid intrinsics
*/
std::string printOverlapInformation(const OverlapMap& sfm_filenames_indexes){
    std::string result = "Name, Size, New Image, Indexes\n";

    for(const auto & p: sfm_filenames_indexes){
        const std::pair<bool,std::vector<IndexT>> & info = p.second;

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
    OverlapMap& sfm_filenames_indexes){
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

      auto it = sfm_filenames_indexes.find(filename);
      if(it != sfm_filenames_indexes.end()){
        it->second.second.push_back(index);
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

bool getVecs2Align(const openMVG::sfm::SfM_Data& first,
  const openMVG::sfm::SfM_Data& second,
  std::vector<openMVG::Vec3> *first_vecs,
  std::vector<openMVG::Vec3> *second_vecs,
  OverlapMap& sfm_filenames_indexes)
{
  IndexT counter=0;
  std::vector<Pair> unused_cameras;

  for(const auto & p: sfm_filenames_indexes){
    const std::pair<bool,std::vector<IndexT>> & info = p.second;
    // lets make sure we have index in each SfM scene
    if(info.second.size()<2){continue;}

    IndexT p1 = info.second[0];
    IndexT p2 = info.second[1];

    //OPENMVG_LOG_INFO << p1 << "," << p2;

    const auto it1 = first.views.find(p1);
    const auto it2 = second.views.find(p2);
    if (it1 == first.views.end() || it2 == second.views.end()){
      // view id referenced by the overlap map no longer exists in this scene
      unused_cameras.push_back(std::make_pair(p1,p2));
      continue;
    }

    const View * view1 = it1->second.get();
    const View * view2 = it2->second.get();

    if(!first.IsPoseAndIntrinsicDefined(view1)){
      OPENMVG_LOG_INFO << view1->id_view << " has no pose associated in first";
      unused_cameras.push_back(std::make_pair(p1,p2));
    }
    else if(!second.IsPoseAndIntrinsicDefined(view2)){
      OPENMVG_LOG_INFO << view2->id_view << " has no pose associated in second";
      unused_cameras.push_back(std::make_pair(p1,p2));
    }
    else{
      first_vecs -> push_back(first.poses.at(view1->id_pose).center());
      second_vecs -> push_back( second.poses.at(view2->id_pose).center() );
      counter++;
    }
  }

  if(first_vecs->size() == second_vecs->size() && first_vecs->size()>3){
    return true;
  }

  return false;
}

/// Compute a 7DOF rigid transform between the two camera trajectories
/// Uses an iterative trimmed fit: fit, drop the worst-residual correspondences,
/// refit, so that a handful of mis-registered overlapping views don't skew the
/// whole alignment.
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
  if (vec_camPosGT.size() < 4) {
    OPENMVG_LOG_ERROR << "Not enough correspondences to compute a similarity transform";
    return false;
  }

  std::vector<size_t> inlier_indices(vec_camPosGT.size());
  std::iota(inlier_indices.begin(), inlier_indices.end(), 0);

  double S = 1.0;
  Vec3 t = Vec3::Zero();
  Mat3 R = Mat3::Identity();

  constexpr int kMaxIterations = 5;
  constexpr double kKeepRatio = 0.8; // keep the best 80% of correspondences each round
  constexpr size_t kMinCorrespondences = 4;

  for (int iteration = 0; iteration < kMaxIterations; ++iteration)
  {
    // Move input point in appropriate container
    openMVG::Mat x1(3, inlier_indices.size());
    openMVG::Mat x2(3, inlier_indices.size());
    for (size_t i = 0; i < inlier_indices.size(); ++i) {
      x1.col(i) = vec_camPosComputed[inlier_indices[i]];
      x2.col(i) = vec_camPosGT[inlier_indices[i]];
    }
    // Compute rigid transformation p'i = S R pi + t
    openMVG::geometry::FindRTS(x1, x2, &S, &t, &R);
    openMVG::geometry::Refine_RTS(x1, x2, &S, &t, &R);

    if (iteration + 1 == kMaxIterations || inlier_indices.size() <= kMinCorrespondences)
      break;

    // Rank correspondences by how well they agree with the current transform and
    // keep only the best kKeepRatio of them for the next round.
    std::vector<std::pair<double,size_t>> residuals;
    residuals.reserve(inlier_indices.size());
    for (const size_t idx : inlier_indices) {
      const Vec3 predicted = S * R * vec_camPosComputed[idx] + t;
      residuals.emplace_back((predicted - vec_camPosGT[idx]).norm(), idx);
    }
    std::sort(residuals.begin(), residuals.end());

    const size_t keep_count = std::max(kMinCorrespondences,
      static_cast<size_t>(residuals.size() * kKeepRatio));
    if (keep_count >= inlier_indices.size()) break; // nothing left to trim

    inlier_indices.clear();
    inlier_indices.reserve(keep_count);
    for (size_t i = 0; i < keep_count; ++i) {
      inlier_indices.push_back(residuals[i].second);
    }
  }

  vec_camPosComputed_T.resize(vec_camPosGT.size());
  for (size_t i = 0; i  < vec_camPosGT.size(); ++i)
  {
    vec_camPosComputed_T[i] = S * R * ( vec_camPosComputed[i]) + t;
  }

  OPENMVG_LOG_INFO << "Similarity transform used " << inlier_indices.size() << "/"
    << vec_camPosGT.size() << " correspondences after trimming";

  *Sout = S;
  *Rout = R;
  *tout = t;
  return true;
}

bool mergeSfMScenes(openMVG::sfm::SfM_Data& sfm_data, openMVG::sfm::SfM_Data& second_sfm_data,
   const double S, const openMVG::Mat3 R, const openMVG::Vec3 T,
   const OverlapMap& sfm_filenames_indexes,
   const std::unordered_map<IndexT, IndexT>& intrinsic_id_remap,
   double landmark_fusion_pixel_tolerance)
{
    std::set<Pair> common_ids;
    std::set<openMVG::IndexT> second_overlap_ids;
    for(const auto & pair: sfm_filenames_indexes){
        const std::pair<bool,std::vector<IndexT>> & info = pair.second;
        // lets make sure we have index in each SfM scene
        if(info.second.size()<2){continue;}

        IndexT p1 = info.second[0];
        IndexT p2 = info.second[1];

        common_ids.insert( Pair(p1,p2) );
        second_overlap_ids.insert(p2);
    }

    // O(1) remaps built once, instead of scanning common_ids / new_view_pairings
    // per-observation (which is what made merging large scenes slow).
    std::unordered_map<IndexT, IndexT> second_to_first_overlap; // second view id -> first view id
    second_to_first_overlap.reserve(common_ids.size());
    for (const auto & p : common_ids){
        second_to_first_overlap[p.second] = p.first;
    }

    // references
    Views & first_views = sfm_data.views;
    Poses & first_poses = sfm_data.poses;

    Views & second_views = second_sfm_data.views;
    Poses & second_poses = second_sfm_data.poses;

    IndexT first_views_size = first_views.size();

    std::set<IndexT> remove_track_ids = std::set<IndexT>();

    std::unordered_map<IndexT, IndexT> new_view_pairings; // second view id -> new first view id

    // This is where we add the poses from the second scene

    for(auto & iterV : second_views){
        ViewPriors *prior = dynamic_cast<sfm::ViewPriors*>(iterV.second.get());
        IndexT id_view = prior->id_view;
        IndexT id_pose = prior->id_pose;

        // if the pose for the first scene exists then we skip, if the pose doesn't we use
        // the second scenes pose and update the firsts
        if(second_overlap_ids.find(id_view) != second_overlap_ids.end()){
          const IndexT first_id = second_to_first_overlap.at(id_view);
          const View * view1 = sfm_data.views.at(first_id).get();
          const View * view2 = second_sfm_data.views.at(id_view).get();

          if(!sfm_data.IsPoseAndIntrinsicDefined(view1) && second_sfm_data.IsPoseAndIntrinsicDefined(view2)){
              OPENMVG_LOG_INFO << "Pose reinstiated from second sfm scene " ;

              Pose3 pose = second_poses.at(view2->id_pose);
              Vec3 nloc = S * R * ( pose.center() ) + T; // update the camera position to the reference scene
              Pose3 npose = Pose3(pose.rotation(),nloc);

              first_poses[view1->id_pose] = npose;
          }
          continue;
        }

        // need to store the new view id to modify the observation ids
        new_view_pairings[id_view] = first_views_size;

        // if there's a pose then modify the pose position
        if (second_sfm_data.IsPoseAndIntrinsicDefined(prior)){
            // image was not used in reconstruction
            // use the similarity transform and generate the new poses
            Pose3 pose = second_poses.at(id_pose);
            Vec3 nloc = S * R * ( pose.center() ) + T; // update the camera position to the reference scene
            Pose3 npose = Pose3(pose.rotation(),nloc);

            first_poses[first_views_size] = npose;

        }else{
            remove_track_ids.insert(id_view);
        }

        prior->id_view = first_views_size;// new view id
        prior->id_pose = first_views_size;// new pose id
        // remap the intrinsic id into the merged scene's intrinsic id-space so this
        // view keeps pointing at the correct (possibly renumbered) camera
        const auto intrinsic_it = intrinsic_id_remap.find(prior->id_intrinsic);
        if (intrinsic_it != intrinsic_id_remap.end()){
            prior->id_intrinsic = intrinsic_it->second;
        }

        first_views[first_views_size] = std::make_shared<ViewPriors>(*prior);

        first_views_size++;
    }

    OPENMVG_LOG_INFO << "Views with no poses: " << remove_track_ids.size();

    // Index the merged scene's existing (first scene) landmarks that are observed
    // by an overlap view, keyed by that view id, so incoming second-scene tracks
    // that observe the same physical image can be fused into the existing landmark
    // instead of being inserted as duplicate 3D points.
    std::unordered_set<IndexT> overlap_first_view_ids;
    overlap_first_view_ids.reserve(common_ids.size());
    for (const auto & p : common_ids){
        overlap_first_view_ids.insert(p.first);
    }

    std::unordered_map<IndexT, std::vector<std::pair<Vec2, IndexT>>> overlap_view_observations;
    if (!overlap_first_view_ids.empty()){
        for (const auto & track : sfm_data.GetLandmarks()){
            for (const auto & obs : track.second.obs){
                if (overlap_first_view_ids.find(obs.first) != overlap_first_view_ids.end()){
                    overlap_view_observations[obs.first].emplace_back(obs.second.x, track.first);
                }
            }
        }
    }
    const double sq_tolerance = landmark_fusion_pixel_tolerance * landmark_fusion_pixel_tolerance;

    const Landmarks & second_tracks = second_sfm_data.GetLandmarks();

    IndexT new_track_counter = sfm_data.structure.size();
    IndexT original_track_num = sfm_data.structure.size();
    IndexT fused_track_count = 0;

    // first update the landmarks in the second scene to the reference frame of the first
    for (const auto& track: second_tracks)
    {
        const Landmark & landmark = track.second;

        Observations new_observations;
        IndexT fusion_target = UndefinedIndexT;

        for (const auto& iterOb: landmark.obs)
        {
            IndexT id_view = iterOb.first;
            // need to update the view_id to the most up to date
            if(second_overlap_ids.find(id_view) != second_overlap_ids.end()){
                // check if the second view has any new information
                const auto it = second_to_first_overlap.find(id_view);
                if (it == second_to_first_overlap.end()){ continue; }
                const IndexT first_view_id = it->second;

                // does an existing landmark already observe this same physical image at
                // (approximately) the same pixel location? if so this is the same 3D point.
                if (fusion_target == UndefinedIndexT){
                    const auto view_obs_it = overlap_view_observations.find(first_view_id);
                    if (view_obs_it != overlap_view_observations.end()){
                        for (const auto & candidate : view_obs_it->second){
                            if ((candidate.first - iterOb.second.x).squaredNorm() <= sq_tolerance){
                                fusion_target = candidate.second;
                                break;
                            }
                        }
                    }
                }

                id_view = first_view_id;
            }else{
                // else we want to use the first views pose if present
                const auto it = new_view_pairings.find(id_view);
                if (it == new_view_pairings.end()){ continue; }
                id_view = it->second;
            }

            const auto view_it = sfm_data.views.find(id_view);
            if (view_it == sfm_data.views.end()){
                OPENMVG_LOG_WARNING << "View id does not exist " << id_view;
                continue;
            }
            if(!sfm_data.IsPoseAndIntrinsicDefined(view_it->second.get())){
                OPENMVG_LOG_WARNING << "Pose not defined for view " << id_view;
                continue;
            }

            // view id will have been updated by now, just have to insert it
            // iterOb.second.id_feat
            new_observations.insert({id_view, Observation(iterOb.second.x, UndefinedIndexT)});
            // new observations will have been added
        }

        if (new_observations.empty()){
            continue;
        }

        const Vec3 new_position = S * R * ( landmark.X ) + T;

        if (fusion_target != UndefinedIndexT){
            // fuse into the existing landmark rather than creating a duplicate 3D point:
            // keep the existing point position, only add observations for views it
            // didn't already have.
            Landmark & target = sfm_data.structure.at(fusion_target);
            for (const auto & obs : new_observations){
                target.obs.emplace(obs.first, obs.second);
            }
            ++fused_track_count;
            continue;
        }

        Landmark new_landmark;
        new_landmark.X = new_position;
        new_landmark.obs = std::move(new_observations);

        sfm_data.structure[new_track_counter++] = std::move( new_landmark );
    }

    OPENMVG_LOG_INFO << "Landmarks fused into existing overlap points: " << fused_track_count
      << ", new landmarks added: " << (new_track_counter - original_track_num);

    return true;
}

}//end of namespace sfm
}//end of namespace openmvg
