/* 
 * File:   optimization.cpp
 * Author: sgaspari
 * 
 * Created on October 23, 2015, 12:02 AM
 */

#include "optimization.hpp"
#include <openMVG/sfm/sfm_data_BA_ceres.hpp>

//@fixme move/redefine
#define POPART_COUT(x) std::cout << x << std::endl
#define POPART_CERR(x) std::cerr << x << std::endl

namespace openMVG{
namespace localization{

bool refineSequence(cameras::Pinhole_Intrinsic_Radial_K3 *intrinsics,
                    std::vector<LocalizationResult> & localizationResult)
{
  
  // flags for BA
  //@todo expose as parameters?
  const bool b_refine_pose = true;
  const bool b_refine_intrinsic = true;
  const bool b_refine_structure = false;
  
  const size_t numCameras = localizationResult.size();
  
  // the id for the instrinsic group
  const IndexT intrinsicID = 0;
    
  // Setup a tiny SfM scene with the corresponding 2D-3D data
  sfm::SfM_Data tinyScene;
  
  for(size_t viewID = 0; viewID < numCameras; ++viewID)
  {
    const LocalizationResult &currResult = localizationResult[viewID];
    // skip invalid poses
    if(!currResult.isValid())
    {
      std::cout << "\n*****\nskipping invalid View " << viewID << std::endl;
      continue;
    }
    
    std::cout << "\n*****\nView " << viewID << std::endl;
    // view
    tinyScene.views.insert( std::make_pair(viewID, std::make_shared<sfm::View>("",viewID, intrinsicID, viewID)));
    // pose
    tinyScene.poses[viewID] = currResult.getPose();
    // intrinsic (the shared_ptr does not take the ownership, will not release the input pointer)
    tinyScene.intrinsics[intrinsicID] = std::shared_ptr<cameras::Pinhole_Intrinsic_Radial_K3>(intrinsics, [](cameras::Pinhole_Intrinsic_Radial_K3*){});
    
    // structure data (2D-3D correspondences)
    const vector<pair<IndexT, IndexT> > &currentIDs = currResult.getIndMatch3D2D();
    
    for(const size_t idx : currResult.getInliers() )
    {
      // the idx should be in the size range of the data
      assert(idx < currResult.getPt3D().cols());
      // get the corresponding 3D point (landmark) ID
      const IndexT landmarkID = currentIDs[idx].first;
      // get the corresponding 2D point ID
      const IndexT featID = currentIDs[idx].second;
      std::cout << "inlier " << idx << " is land " << landmarkID << " and feat " << featID << std::endl;
      // get the corresponding feature
      const Vec2 &feature = currResult.getPt2D().col(idx);
      // check if the point exists already
      if(tinyScene.structure.count(landmarkID))
      {
        // normally there should be no other features already associated to this
        // 3D point in this view
//        assert(tinyScene.structure[landmarkID].obs.count(viewID) == 0);
        if(tinyScene.structure[landmarkID].obs.count(viewID) != 0)
        {
          // this is weird but it could happen when two features are really close to each other (?)
          std::cout << "Point 3D " << landmarkID << " has multiple features " << tinyScene.structure[landmarkID].obs.size() << " in the same view " << viewID << " size "  << std::endl; 
          continue;
        }
        
        // the 3D point exists already, add the observation
        tinyScene.structure[landmarkID].obs[viewID] =  sfm::Observation(feature, featID);
      }
      else
      {
        // create a new 3D point
        sfm::Landmark landmark;
        landmark.X = currResult.getPt3D().col(idx);
        landmark.obs[viewID] = sfm::Observation(feature, featID);
        tinyScene.structure[landmarkID] = std::move(landmark);
      }
    }
  }
  POPART_COUT("Number of 3D-2D associations " << tinyScene.structure.size());
  std::vector<double> params = intrinsics->getParams();
  POPART_COUT("K before bundle:" << params[0] << " " << params[1] << " "<< params[2]);
  POPART_COUT("Distortion before bundle" << params[3] << " " << params[4] << " "<< params[5]);

  sfm::Bundle_Adjustment_Ceres bundle_adjustment_obj;
  const bool b_BA_Status = bundle_adjustment_obj.Adjust(tinyScene, b_refine_pose, b_refine_pose, b_refine_intrinsic, b_refine_structure);
  if(b_BA_Status)
  {
    // get back the results
    for(const auto &pose : tinyScene.poses)
    {
      const IndexT idPose = pose.first;
      localizationResult[idPose].setPose(pose.second);
    }
  }
  params = intrinsics->getParams();
  POPART_COUT("K after bundle:" << params[0] << " " << params[1] << " "<< params[2]);
  POPART_COUT("Distortion after bundle" << params[3] << " " << params[4] << " "<< params[5]);
  
  return b_BA_Status;
}

} //namespace localization
} //namespace openMVG

