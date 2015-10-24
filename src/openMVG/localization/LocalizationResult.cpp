#include "LocalizationResult.hpp"

namespace openMVG {
namespace localization {

LocalizationResult::LocalizationResult() : 
        _isValid(false)
{
}

LocalizationResult::LocalizationResult(
        const sfm::Image_Localizer_Match_Data & matchData,
        const std::vector<pair<IndexT, IndexT> > & indMatch3D2D,
        const geometry::Pose3 & pose,
        const cameras::Pinhole_Intrinsic_Radial_K3 & intrinsics,
        bool isValid) :
        _matchData(matchData),
        _indMatch3D2D(indMatch3D2D),
        _pose(pose),
        _intrinsics(intrinsics),
        _isValid(isValid)
{
}

LocalizationResult::~LocalizationResult()
{
}

// Accessor
const std::vector<size_t> & LocalizationResult::getInliers() const 
{
  return _matchData.vec_inliers;
}

const Mat & LocalizationResult::getPt2D() const 
{
  return _matchData.pt2D;
}

const Mat & LocalizationResult::getPt3D() const 
{
  return _matchData.pt3D;
}

const std::vector<pair<IndexT, IndexT> > & LocalizationResult::getIndMatch3D2D() const
{
  return _indMatch3D2D;
}
        
const geometry::Pose3 & LocalizationResult::getPose() const
{
  return _pose;
}

void LocalizationResult::setPose(const geometry::Pose3 & pose)
{
  _pose = pose;
}

const cameras::Pinhole_Intrinsic_Radial_K3 & LocalizationResult::getIntrinsics() const
{
  return _intrinsics;
}

bool LocalizationResult::isValid() const
{
  return _isValid;
}

} // localization
} // openMVG