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
        const std::vector<voctree::DocMatch>& matchedImages,
        bool isValid) :
        _matchData(matchData),
        _indMatch3D2D(indMatch3D2D),
        _pose(pose),
        _intrinsics(intrinsics),
        _matchedImages(matchedImages),
        _isValid(isValid)
{
  // verify data consistency
  assert(_matchData.pt2D.cols() == _matchData.pt3D.cols());
  assert(_matchData.pt2D.cols() == _indMatch3D2D.size());
}
        
LocalizationResult::~LocalizationResult()
{
}

const Mat LocalizationResult::retrieveUndistortedPt2D() const
{
  const auto &intrinsics =  getIntrinsics();
  const auto &distorted = getPt2D();
  if(!intrinsics.have_disto() || !intrinsics.isValid())
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

Mat2X LocalizationResult::computeAllResiduals() const 
{
  const Mat2X &orig2d = getPt2D();
  const Mat3X &orig3d = getPt3D();
  assert(orig2d.cols()==orig3d.cols());
  
  const auto &intrinsics = getIntrinsics();
  return intrinsics.residuals(getPose(), orig3d, orig2d);
}

Mat2X LocalizationResult::computeInliersResiduals() const 
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

double LocalizationResult::computeInliersRMSE() const 
{
  const auto& residuals = computeInliersResiduals();
  // squared residual for each point
  const auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
  //RMSE
  return std::sqrt(sqrErrors.mean());
}

double LocalizationResult::computeAllRMSE() const 
{
  const auto& residuals = computeAllResiduals();
  // squared residual for each point
  const auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
  //RMSE
  return std::sqrt(sqrErrors.mean());
}

Vec LocalizationResult::computeReprojectionErrorPerInlier() const 
{
  const auto& residuals = computeInliersResiduals();
  // squared residual for each point
  const Vec sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
  //RMSE
  return sqrErrors.array().sqrt();
}

Vec LocalizationResult::computeReprojectionErrorPerPoint() const
{
  const auto& residuals = computeAllResiduals();
  // squared residual for each point
  const Vec sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
  //RMSE
  return sqrErrors.array().sqrt();
}

std::size_t LocalizationResult::selectBestInliers(double maxReprojectionError)
{
   const auto &residuals = computeReprojectionErrorPerPoint();
   auto &inliers = _matchData.vec_inliers;
   std::cout << "Inliers before: " << inliers.size();
   inliers.clear();
   // at worst they could all be inliers
   inliers.reserve(getPt2D().size());
   
   for(std::size_t i = 0; i < residuals.size(); ++i )
   {
     if(residuals[i] < maxReprojectionError)
     {
       inliers.push_back(i);
     }
   }
   std::cout << " \tafter: " << inliers.size() << std::endl;
   _matchData.error_max = maxReprojectionError;
   return inliers.size();
}

std::size_t LocalizationResult::selectBestInliers()
{  
   const auto threshold = _matchData.error_max;
   return selectBestInliers(threshold);
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



void updateRigPoses(std::vector<LocalizationResult>& vec_localizationResults,
                    const geometry::Pose3 &rigPose,
                    const std::vector<geometry::Pose3 > &vec_subPoses)
{
  const std::size_t numCams = vec_localizationResults.size();
  assert(numCams==vec_subPoses.size()+1);
  
  // update localization result poses
  for(std::size_t camID = 0; camID < numCams; ++camID)
  {
    geometry::Pose3 pose;
    if(camID == 0)
    {
      pose = rigPose;
    }
    else
    {
      // main camera: q1 ~ [R1 t1] Q = [I 0] A   where A = [R1 t1] Q  
      // another camera: q2 ~ [R2 t2] Q = [R2 t2]*inv([R1 t1]) A 
      // and subPose12 = [R12 t12] = [R2 t2]*inv([R1 t1])
      // With rigResection() we compute [R1 t1] (aka rigPose), hence:
      // subPose12 = [R12 t12] = [R2 t2]*inv([R1 t1]) and we need [R2 t2], ie the absolute pose
      // => [R1 t1] * subPose12 = [R2 t2]
      // => rigPose * subPose12 = [R2 t2]
      pose = vec_subPoses[camID-1] * rigPose;
    }
    
    vec_localizationResults[camID].setPose(pose);
  }
}

} // localization
} // openMVG