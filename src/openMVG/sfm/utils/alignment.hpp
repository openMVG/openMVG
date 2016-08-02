#pragma once

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/geometry/pose3.hpp"

namespace openMVG {
namespace sfm {

inline void getCommonViews(const SfM_Data & sfmDataA,
                           const SfM_Data & sfmDataB,
                           std::vector<IndexT>& outIndexes)
{
  for(const auto& viewA: sfmDataA.GetViews())
  {
    if(sfmDataB.GetViews().find(viewA.first) != sfmDataB.GetViews().end())
    {
      outIndexes.push_back(viewA.first);
    }
  }
}

inline void getCommonViewsWithPoses(const SfM_Data & sfmDataA,
                                    const SfM_Data & sfmDataB,
                                    std::vector<IndexT>& outIndexes)
{
  for(const auto& viewA: sfmDataA.GetViews())
  {
    // check there is a view with the same ID and both of them have pose and 
    // intrinsics defined
    if(!sfmDataA.IsPoseAndIntrinsicDefined(viewA.second.get()))
      continue;

    if(sfmDataB.GetViews().find(viewA.first) != sfmDataB.GetViews().end() &&
       sfmDataB.IsPoseAndIntrinsicDefined(viewA.first))
    {
      outIndexes.push_back(viewA.first);
    }
  }
}

/**
 * @brief Compute a 5DOF rigid transform between the two set of cameras.
 *
 * @param[in] sfmDataA
 * @param[in] sfmDataB
 * @param[out] out_S output scale factor
 * @param[out] out_R output rotation 3x3 matrix
 * @param[out] out_t output translation vector
 * @return true if it finds a similarity transformation
 */
bool computeSimilarity(const SfM_Data & sfmDataA,
                       const SfM_Data & sfmDataB,
                       double * out_S,
                       Mat3 * out_R,
                       Vec3 * out_t);


inline void applyTransform(SfM_Data & sfmData,
                           const double S,
                           const Mat3& R,
                           const Vec3& t)
{
  for(auto& view: sfmData.views)
  {
    geometry::Pose3& pose = sfmData.poses[view.second->id_pose];
    pose = pose.transformSRt(S, R, t);
  }
  for(auto& landmark: sfmData.structure)
  {
    landmark.second.X = S * R * landmark.second.X + t;
  }
}

}
}
