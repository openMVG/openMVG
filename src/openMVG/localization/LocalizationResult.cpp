#include "LocalizationResult.hpp"

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/utility.hpp>  // needed to serialize std::pair

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

const Mat LocalizationResult::getUndistortedPt2D() const
{
  const auto &intrinsics =  getIntrinsics();
  const auto &distorted = getPt2D();
  if(!intrinsics.have_disto())
  {
    return getPt2D();
  }
  const std::size_t numPts = distorted.cols();
  Mat pt2Dundistorted = Mat2X(2, numPts);
  for(std::size_t iPoint = 0; iPoint < numPts; ++iPoint)
  {
    pt2Dundistorted.col(iPoint) = intrinsics.get_ud_pixel(distorted.col(iPoint));
  }
  return pt2Dundistorted;
}

const Mat & LocalizationResult::getPt3D() const 
{
  return _matchData.pt3D;
}

 const Mat34 & LocalizationResult::getProjection() const
 {
   return _matchData.projection_matrix;
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

cameras::Pinhole_Intrinsic_Radial_K3 & LocalizationResult::getIntrinsics()
{
  return _intrinsics;
}

void LocalizationResult::updateIntrinsics(const std::vector<double> & params)
{
  _intrinsics.updateFromParams(params);
}

bool LocalizationResult::isValid() const
{
  return _isValid;
}

Mat2X LocalizationResult::computeResiduals() const 
{
  // get the inliers.
  const auto &currInliers = getInliers();
  const std::size_t numInliers = currInliers.size();

  const Mat2X &orig2d = getPt2D();
  const Mat3X &orig3d = getPt3D();
  Mat2X inliers2d = Mat2X(2, numInliers);
  Mat3X inliers3d = Mat3X(3, numInliers);

  for(std::size_t i = 0; i < numInliers; ++i)
  {
    const std::size_t idx = currInliers[i];
    inliers2d.col(i) = orig2d.col(idx);
    inliers3d.col(i) = orig3d.col(idx);
  }
  
  const auto &intrinsics = getIntrinsics();
  return intrinsics.residuals(getPose(), inliers3d, inliers2d);
}



bool load(LocalizationResult & res, const std::string & filename)
{
  //Create the stream and check it is ok
  std::ifstream stream(filename, std::ios::binary | std::ios::in);
  if(!stream.is_open())
  {
    std::cerr << "Unable to load file " << filename << std::endl;
    return false;
  }
  try
  {
    cereal::PortableBinaryInputArchive archive(stream);
    archive(cereal::make_nvp("result", res));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  return true;
}


bool load(std::vector<LocalizationResult> & res, const std::string & filename)
{
  //Create the stream and check it is ok
  std::ifstream stream(filename, std::ios::binary | std::ios::in);
  if(!stream.is_open())
  {
    std::cerr << "Unable to load file " << filename << std::endl;
    return false;
  }
  try
  {
    cereal::PortableBinaryInputArchive archive(stream);
    archive(cereal::make_nvp("results", res));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  return true;
}


bool save(const LocalizationResult & res, const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename, std::ios::binary | std::ios::out);
  if(!stream.is_open())
  {
    std::cerr << "Unable to create file " << filename << std::endl;
    return false;
  }

  cereal::PortableBinaryOutputArchive archive(stream);
  archive(res);

  return true; 
}

bool save(const std::vector<LocalizationResult> & res, const std::string & filename)
{
  //Create the stream and check it is ok
  std::ofstream stream(filename, std::ios::binary | std::ios::out);
  if(!stream.is_open())
  {
    std::cerr << "Unable to create file " << filename << std::endl;
    return false;
  }

  cereal::PortableBinaryOutputArchive archive(stream);
  archive(cereal::make_nvp("results", res));

  return true;  
}


} // localization
} // openMVG