#pragma once

#include <openMVG/localization/localization.hpp>

#include "openMVG/features/features.hpp"
#include <openMVG/features/image_describer.hpp>
#include <openMVG/features/cctag/CCTAG_describer.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/pipelines/localization/SfM_Localizer.hpp>

#include <iostream>

namespace openMVG {
namespace localization {
  
typedef Reconstructed_RegionsT::DescriptorT DescriptorT;

class CCTagLocalizer {
public:
  
  bool init(const std::string &sfmFilePath,
            const std::string &descriptorsFolder);
  
  bool Localize(const image::Image<unsigned char> & imageGray,
                cameras::IntrinsicBase * queryIntrinsics,
                geometry::Pose3 & pose,
                bool useGuidedMatching,
                sfm::Image_Localizer_Match_Data * resection_data = nullptr);
  
  CCTagLocalizer();
  CCTagLocalizer(const CCTagLocalizer& orig);
  virtual ~CCTagLocalizer();
private:
  
  bool loadReconstructionDescriptors(
    const sfm::SfM_Data & sfm_data,
    const std::string & feat_directory);
  
  // for each view index, it contains the cctag features and descriptors that have an
  // associated 3D point
  Hash_Map<IndexT, Reconstructed_RegionsT > _regions_per_view;
  
  // contains the 3D reconstruction data
  sfm::SfM_Data _sfm_data;
  
  // the feature extractor
  features::CCTAG_Image_describer _image_describer;
  
  //
  std::map<IndexT, Vec3> _cctagDatabase;
};

IndexT getCCTagId(const DescriptorT & desc);

} // namespace localization
} // openMVG
