#include "LocalizationResult.hpp"
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>
#include <openMVG/numeric/numeric.h>
#include <openMVG/sfm/pipelines/localization/SfM_Localizer.hpp>
#include "testing/testing.h"

#include <third_party/stlplus3/filesystemSimplified/file_system.hpp>

#include <vector>
#include <chrono>
#include <random>

using namespace openMVG;

sfm::Image_Localizer_Match_Data generateRandomMatch_Data(std::size_t numPts)
{
  sfm::Image_Localizer_Match_Data data;
  data.projection_matrix = Mat34::Random();
  data.pt3D = Mat::Random(3, numPts);
  data.pt2D = Mat::Random(2, numPts);
  const std::size_t numInliers = (std::size_t) numPts*0.7;
  for(std::size_t i = 0; i < numInliers; ++i)
  {
    data.vec_inliers.push_back(i);
  }
  return data;
}

localization::LocalizationResult generateRandomResult(std::size_t numPts)
{
  // random matchData
  const sfm::Image_Localizer_Match_Data &data = generateRandomMatch_Data(numPts);
  
  // random indMatch3D2D
  std::vector<pair<IndexT, IndexT> > indMatch3D2D;
  indMatch3D2D.reserve(numPts);
  for(std::size_t i = 0; i < numPts; ++i)
  {
    indMatch3D2D.emplace_back(i,i);
  }
  
  // random pose
  geometry::Pose3 pose = geometry::Pose3(Mat3::Random(), Vec3::Random());
  
  // random intrinsics
  cameras::Pinhole_Intrinsic_Radial_K3 intrinsics = cameras::Pinhole_Intrinsic_Radial_K3(640, 480, 1400, 320.5, 240.5, 0.001, -0.05, 0.00003);
  
  // random valid
  const bool valid = (numPts % 2 == 0);

  std::vector<voctree::DocMatch> matchedImages;
  matchedImages.push_back(voctree::DocMatch(2, 0.5));
  matchedImages.push_back(voctree::DocMatch(3, 0.8));

  return localization::LocalizationResult(data, indMatch3D2D, pose, intrinsics, matchedImages, valid);
}

// generate a random localization result, save it to binary file, load it again
// and compare each value

TEST(LocalizationResult, LoadSaveBinSingle)
{
  const double threshold = 1e-10;
  const std::string filename = "test_localizationResult.bin";
  const localization::LocalizationResult &res = generateRandomResult(10);
  localization::LocalizationResult check;

  EXPECT_TRUE(localization::save(res, filename));
  EXPECT_TRUE(localization::load(check, filename));

  // same validity
  EXPECT_TRUE(res.isValid() == check.isValid());

  // same pose
  const Mat3 rotGT = res.getPose().rotation();
  const Mat3 rot = check.getPose().rotation();
  for(std::size_t i = 0; i < 3; ++i)
  {
    for(std::size_t j = 0; j < 3; ++j)
    {
      EXPECT_NEAR(rotGT(i, j), rot(i, j), threshold)
    }
  }
  const Vec3 centerGT = res.getPose().center();
  const Vec3 center = check.getPose().center();
  EXPECT_NEAR(centerGT(0), center(0), threshold);
  EXPECT_NEAR(centerGT(1), center(1), threshold);
  EXPECT_NEAR(centerGT(2), center(2), threshold);

  // same _indMatch3D2D
  const auto idxGT = res.getIndMatch3D2D();
  const auto idx = check.getIndMatch3D2D();
  EXPECT_TRUE(idxGT.size() == idx.size());
  const std::size_t numpts = idxGT.size();
  for(std::size_t i = 0; i < numpts; ++i)
  {
    EXPECT_TRUE(idxGT[i].first == idx[i].first);
    EXPECT_TRUE(idxGT[i].second == idx[i].second);
  }

  // same _matchData
  EXPECT_TRUE(res.getInliers().size() == check.getInliers().size());
  const auto inliersGT = res.getInliers();
  const auto inliers = check.getInliers();
  for(std::size_t i = 0; i < res.getInliers().size(); ++i)
  {
    EXPECT_TRUE(inliersGT[i] == inliers[i]);
  }


  EXPECT_MATRIX_NEAR(res.getPt3D(), check.getPt3D(), threshold);
  EXPECT_MATRIX_NEAR(res.getPt2D(), check.getPt2D(), threshold);
  EXPECT_MATRIX_NEAR(res.getProjection(), check.getProjection(), threshold);

    // same matchedImages
  EXPECT_TRUE(res.getMatchedImages().size() == check.getMatchedImages().size());
  const std::vector<voctree::DocMatch>& matchedImagesGT = res.getMatchedImages();
  const std::vector<voctree::DocMatch>& matchedImages = check.getMatchedImages();
  for(std::size_t i = 0; i < res.getMatchedImages().size(); ++i)
  {
    EXPECT_TRUE(matchedImagesGT[i] == matchedImages[i]);
  }

  stlplus::file_delete(filename);
}

TEST(LocalizationResult, LoadSaveBinVector)
{
  const double threshold = 1e-10;
  const std::size_t numResults = 10;
  const std::string filename = "test_localizationResults.bin";
  const unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  std::random_device rd;
  std::mt19937 gen(seed1);
  std::uniform_int_distribution<> numpts(1, 20);

  std::vector<localization::LocalizationResult> resGT;
  std::vector<localization::LocalizationResult> resCheck;
  resGT.reserve(numResults);
  resCheck.reserve(numResults);

  for(std::size_t i = 0; i < numResults; ++i)
  {
    resGT.push_back(generateRandomResult(numpts(gen)));
  }

  EXPECT_TRUE(localization::save(resGT, filename));
  EXPECT_TRUE(localization::load(resCheck, filename));
  EXPECT_TRUE(resCheck.size() == resGT.size());

  // check each element
  for(std::size_t i = 0; i < numResults; ++i)
  {
    const auto res = resGT[i];
    const auto check = resCheck[i];

    // same validity
    EXPECT_TRUE(res.isValid() == check.isValid());

    // same pose
    const Mat3 rotGT = res.getPose().rotation();
    const Mat3 rot = check.getPose().rotation();
    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
      {
        EXPECT_NEAR(rotGT(i, j), rot(i, j), threshold)
      }
    }
    const Vec3 centerGT = res.getPose().center();
    const Vec3 center = check.getPose().center();
    EXPECT_NEAR(centerGT(0), center(0), threshold);
    EXPECT_NEAR(centerGT(1), center(1), threshold);
    EXPECT_NEAR(centerGT(2), center(2), threshold);

    // same _indMatch3D2D
    const auto idxGT = res.getIndMatch3D2D();
    const auto idx = check.getIndMatch3D2D();
    EXPECT_TRUE(idxGT.size() == idx.size());
    const std::size_t numpts = idxGT.size();
    for(std::size_t j = 0; j < numpts; ++j)
    {
      EXPECT_TRUE(idxGT[j].first == idx[j].first);
      EXPECT_TRUE(idxGT[j].second == idx[j].second);
    }

    // same _matchData
    EXPECT_TRUE(res.getInliers().size() == check.getInliers().size());
    const auto inliersGT = res.getInliers();
    const auto inliers = check.getInliers();
    for(std::size_t j = 0; j < res.getInliers().size(); ++j)
    {
      EXPECT_TRUE(inliersGT[j] == inliers[j]);
    }

    EXPECT_MATRIX_NEAR(res.getPt3D(), check.getPt3D(), threshold);
    EXPECT_MATRIX_NEAR(res.getPt2D(), check.getPt2D(), threshold);
    EXPECT_MATRIX_NEAR(res.getProjection(), check.getProjection(), threshold);
    
    // same matchedImages
    EXPECT_TRUE(res.getMatchedImages().size() == check.getMatchedImages().size());
    const auto matchedImagesGT = res.getMatchedImages();
    const auto matchedImages = check.getMatchedImages();
    for(std::size_t j = 0; j < res.getMatchedImages().size(); ++j)
    {
      EXPECT_TRUE(matchedImagesGT[j] == matchedImages[j]);
    }

    stlplus::file_delete(filename);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
