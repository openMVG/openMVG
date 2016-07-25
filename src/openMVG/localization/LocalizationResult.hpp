#pragma once

#include <openMVG/sfm/pipelines/localization/SfM_Localizer.hpp>
#include <cereal/cereal.hpp>
#include <vector>

namespace cereal{

template<class Archive>
void save(Archive & archive, openMVG::Mat34 const & mat)
{ 
  vector<double> vec(mat.data(), mat.data() + mat.rows() * mat.cols());
  archive(vec); 
}

template<class Archive>
void load(Archive & archive, openMVG::Mat34 & m)
{
  vector<double> vec(12);
  archive(vec); 
  m = Eigen::Map<openMVG::Mat34>(vec.data(), 3, 4);
}

}

namespace openMVG {
namespace localization {

template<class Archive>
void saveMat(Archive & archive, const std::string &name, Mat const & mat)
{ 
  vector<double> vec(mat.data(), mat.data() + mat.rows() * mat.cols());
  archive(cereal::make_nvp(name, vec)); 
}

template<class Archive>
void loadMat(Archive & archive, const std::string &name, const std::size_t rows, openMVG::Mat & m)
{
  vector<double> vec;
  archive(cereal::make_nvp(name, vec));
  const size_t cols = vec.size() / rows;
  m = Eigen::Map<openMVG::Mat>(vec.data(), rows, cols);
}

class LocalizationResult {
public:
  
  LocalizationResult();
  
  LocalizationResult(
        const sfm::Image_Localizer_Match_Data & matchData,
        const std::vector<pair<IndexT, IndexT> > & indMatch3D2D,
        const geometry::Pose3 & pose,
        const cameras::Pinhole_Intrinsic_Radial_K3 & intrinsics,
        bool isValid = true);
  
  virtual ~LocalizationResult();
  
  // Accessors
  const std::vector<size_t> & getInliers() const;

  const Mat & getPt2D() const;
  
  const Mat getUndistortedPt2D() const;

  const Mat & getPt3D() const;
  
  const Mat34 & getProjection() const;
  
  sfm::Image_Localizer_Match_Data getMatchData() const { return _matchData; } 

  const std::vector<pair<IndexT, IndexT> > & getIndMatch3D2D() const;

  const geometry::Pose3 & getPose() const;
  
  void setPose(const geometry::Pose3 & pose);
  
  const cameras::Pinhole_Intrinsic_Radial_K3 & getIntrinsics() const;
  
  cameras::Pinhole_Intrinsic_Radial_K3 & getIntrinsics();
  
  void updateIntrinsics(const std::vector<double> & params);

  bool isValid() const;
  
  Mat2X computeResiduals() const ;

  // Serialization
  template<class Archive>
  void save(Archive & ar) const
  {
    ar(cereal::make_nvp("isValid", _isValid));
    ar(cereal::make_nvp("pose", _pose));
    ar(cereal::make_nvp("indMatch3D2D", _indMatch3D2D));
    ar(cereal::make_nvp("intrinsics", _intrinsics));
    ar(cereal::make_nvp("inliers", _matchData.vec_inliers));
    saveMat(ar, "pt3D", _matchData.pt3D);
    saveMat(ar, "pt2D", _matchData.pt2D);
    ar(cereal::make_nvp("projection_matrix", _matchData.projection_matrix));
    ar(cereal::make_nvp("error_max", _matchData.error_max));
    ar(cereal::make_nvp("max_iteration", _matchData.max_iteration));
  }

  template<class Archive>
  void load(Archive & ar)
  {
    ar(cereal::make_nvp("isValid", _isValid));
    ar(cereal::make_nvp("pose", _pose));
    ar(cereal::make_nvp("indMatch3D2D", _indMatch3D2D));
    ar(cereal::make_nvp("intrinsics", _intrinsics));
    ar(cereal::make_nvp("inliers", _matchData.vec_inliers));
    loadMat(ar, "pt3D", 3, _matchData.pt3D);
    loadMat(ar, "pt2D", 2 , _matchData.pt2D);
    ar(cereal::make_nvp("projection_matrix", _matchData.projection_matrix));
    ar(cereal::make_nvp("error_max", _matchData.error_max));
    ar(cereal::make_nvp("max_iteration", _matchData.max_iteration));
  }
  
private:
  sfm::Image_Localizer_Match_Data _matchData;
                             // Hold all the imaged points, 
                             // their associated 3D points and the inlier indices 
                             // (w.r.t. the F/E robust estimation)

  std::vector<std::pair<IndexT, IndexT> > _indMatch3D2D; // 3D to 2D index matches in the global index system,
                                                    // i.e. the set of pair (landmark id, index of the associated 2D point)

  geometry::Pose3 _pose; // Computed camera pose
  
  cameras::Pinhole_Intrinsic_Radial_K3 _intrinsics;

  bool _isValid; // True if the localization succeeded, false otherwise
  
private:

};

bool load(LocalizationResult & res, const std::string & filename);
bool load(std::vector<LocalizationResult> & res, const std::string & filename);


bool save(const LocalizationResult & res, const std::string & filename);
bool save(const std::vector<LocalizationResult> & res, const std::string & filename);

} // localization
} // openMVG

